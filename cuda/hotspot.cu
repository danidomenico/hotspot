%%writefile hotspot.cu

#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <cassert>
#include <time.h>
#include <sys/time.h>

using FLOAT = float;

// Returns the current system time in microseconds 
long long get_time() {
	struct timeval tv;
	gettimeofday(&tv, NULL);
	return (tv.tv_sec * 1000000) + tv.tv_usec;
}

#ifdef RD_WG_SIZE_0_0                                                            
        #define BLOCK_SIZE RD_WG_SIZE_0_0                                        
#elif defined(RD_WG_SIZE_0)                                                      
        #define BLOCK_SIZE RD_WG_SIZE_0                                          
#elif defined(RD_WG_SIZE)                                                        
        #define BLOCK_SIZE RD_WG_SIZE                                            
#else                                                                                    
        #define BLOCK_SIZE 16                                                            
#endif                                                                                   

#define STR_SIZE 256

/* maximum power density possible (say 300W for a 10mm x 10mm chip)	*/
#define MAX_PD	(3.0e6)
/* required precision in degrees	*/
#define PRECISION	0.001
#define SPEC_HEAT_SI 1.75e6
#define K_SI 100
/* capacitance fitting factor	*/
#define FACTOR_CHIP	0.5

/* chip parameters	*/
float t_chip = 0.0005;
float chip_height = 0.016;
float chip_width = 0.016;
/* ambient temperature, assuming no package at all	*/
float amb_temp = 80.0;

/* define timer macros */
#define pin_stats_reset()   startCycle()
#define pin_stats_pause(cycles)   stopCycle(cycles)
#define pin_stats_dump(cycles)    printf("timer: %Lu\n", cycles)

#define IN_RANGE(x, min, max)   ((x)>=(min) && (x)<=(max))
#define CLAMP_RANGE(x, min, max) x = (x<(min)) ? min : ((x>(max)) ? max : x )
#define MIN(a, b) ((a)<=(b) ? (a) : (b))

__global__ void calculate_temp(int iteration,  //number of iteration
                               float *power,   //power input
                               float *temp_src,    //temperature input/output
                               float *temp_dst,    //temperature input/output
                               int grid_cols,  //Col of grid
                               int grid_rows,  //Row of grid
							   int border_cols,  // border offset 
							   int border_rows,  // border offset
                               float Cap,      //Capacitance
                               float Rx, 
                               float Ry, 
                               float Rz, 
                               float step, 
                               float time_elapsed){
	
        __shared__ float temp_on_cuda[BLOCK_SIZE][BLOCK_SIZE];
        __shared__ float power_on_cuda[BLOCK_SIZE][BLOCK_SIZE];
        __shared__ float temp_t[BLOCK_SIZE][BLOCK_SIZE]; // saving temparary temperature result

	float amb_temp = 80.0;
        float step_div_Cap;
        float Rx_1,Ry_1,Rz_1;
        
	int bx = blockIdx.x;
        int by = blockIdx.y;

	int tx=threadIdx.x;
	int ty=threadIdx.y;
	
	step_div_Cap=step/Cap;
	
	Rx_1=1/Rx;
	Ry_1=1/Ry;
	Rz_1=1/Rz;
	
        // each block finally computes result for a small block
        // after N iterations. 
        // it is the non-overlapping small blocks that cover 
        // all the input data

        // calculate the small block size
	int small_block_rows = BLOCK_SIZE-iteration*2;//EXPAND_RATE
	int small_block_cols = BLOCK_SIZE-iteration*2;//EXPAND_RATE

        // calculate the boundary for the block according to 
        // the boundary of its small block
        int blkY = small_block_rows*by-border_rows;
        int blkX = small_block_cols*bx-border_cols;
        int blkYmax = blkY+BLOCK_SIZE-1;
        int blkXmax = blkX+BLOCK_SIZE-1;

        // calculate the global thread coordination
	int yidx = blkY+ty;
	int xidx = blkX+tx;

        // load data if it is within the valid input range
	int loadYidx=yidx, loadXidx=xidx;
        int index = grid_cols*loadYidx+loadXidx;
       
	if(IN_RANGE(loadYidx, 0, grid_rows-1) && IN_RANGE(loadXidx, 0, grid_cols-1)){
            temp_on_cuda[ty][tx] = temp_src[index];  // Load the temperature data from global memory to shared memory
            power_on_cuda[ty][tx] = power[index];// Load the power data from global memory to shared memory
	}
	__syncthreads();

        // effective range within this block that falls within 
        // the valid range of the input data
        // used to rule out computation outside the boundary.
        int validYmin = (blkY < 0) ? -blkY : 0;
        int validYmax = (blkYmax > grid_rows-1) ? BLOCK_SIZE-1-(blkYmax-grid_rows+1) : BLOCK_SIZE-1;
        int validXmin = (blkX < 0) ? -blkX : 0;
        int validXmax = (blkXmax > grid_cols-1) ? BLOCK_SIZE-1-(blkXmax-grid_cols+1) : BLOCK_SIZE-1;

        int N = ty-1;
        int S = ty+1;
        int W = tx-1;
        int E = tx+1;
        
        N = (N < validYmin) ? validYmin : N;
        S = (S > validYmax) ? validYmax : S;
        W = (W < validXmin) ? validXmin : W;
        E = (E > validXmax) ? validXmax : E;

        bool computed;
        for (int i=0; i<iteration ; i++){ 
            computed = false;
            if( IN_RANGE(tx, i+1, BLOCK_SIZE-i-2) &&  \
                  IN_RANGE(ty, i+1, BLOCK_SIZE-i-2) &&  \
                  IN_RANGE(tx, validXmin, validXmax) && \
                  IN_RANGE(ty, validYmin, validYmax) ) {
                  computed = true;
                  temp_t[ty][tx] =   temp_on_cuda[ty][tx] + step_div_Cap * (power_on_cuda[ty][tx] + 
	       	         (temp_on_cuda[S][tx] + temp_on_cuda[N][tx] - 2.0*temp_on_cuda[ty][tx]) * Ry_1 + 
		             (temp_on_cuda[ty][E] + temp_on_cuda[ty][W] - 2.0*temp_on_cuda[ty][tx]) * Rx_1 + 
		             (amb_temp - temp_on_cuda[ty][tx]) * Rz_1);
	
            }
            __syncthreads();
            if(i==iteration-1)
                break;
            if(computed)	 //Assign the computation range
                temp_on_cuda[ty][tx]= temp_t[ty][tx];
            __syncthreads();
          }

      // update the global memory
      // after the last iteration, only threads coordinated within the 
      // small block perform the calculation and switch on ``computed''
      if (computed){
          temp_dst[index]= temp_t[ty][tx];		
      }
}

/*
   compute N time steps
*/

int compute_tran_temp(FLOAT *MatrixPower, FLOAT *MatrixTemp[2], int col, int row, \
		int total_iterations, int num_iterations, int blockCols, int blockRows, int borderCols, int borderRows) {
    dim3 dimBlock(BLOCK_SIZE, BLOCK_SIZE);
    dim3 dimGrid(blockCols, blockRows);  
	
	float grid_height = chip_height / row;
	float grid_width = chip_width / col;

	float Cap = FACTOR_CHIP * SPEC_HEAT_SI * t_chip * grid_width * grid_height;
	float Rx = grid_width / (2.0 * K_SI * t_chip * grid_height);
	float Ry = grid_height / (2.0 * K_SI * t_chip * grid_width);
	float Rz = t_chip / (K_SI * grid_height * grid_width);

	float max_slope = MAX_PD / (FACTOR_CHIP * t_chip * SPEC_HEAT_SI);
	float step = PRECISION / max_slope;

    float time_elapsed = 0.001;

    int src = 1, dst = 0;
	
	for(int t = 0; t < total_iterations; t+=num_iterations) {
        int temp = src;
        src = dst;
        dst = temp;
        calculate_temp<<<dimGrid, dimBlock>>>(MIN(num_iterations, total_iterations-t), MatrixPower, MatrixTemp[src], MatrixTemp[dst],\
		    col,row, borderCols, borderRows, Cap, Rx, Ry, Rz, step, time_elapsed);
	}
    
    return dst;
}

void fatal(std::string s) {
	std::cerr << "Hotspot error: " << s << std::endl;
	std::abort();
}

void generate_input(FLOAT* temp, FLOAT* power, FLOAT* result, int order_size) {
	const FLOAT MIN = 323.0;
	const FLOAT MAX = 341.0;
	FLOAT value = MIN;
	bool  asc = true;
	
	for(unsigned i=0; i<order_size * order_size; i++) {
		//POWER
		FLOAT x = (rand() % 1000000) / 1000000.0;
		power[i] = x;
		
		//TEMP
		x = (rand() % 1000000) / 1000000.0;
		if(asc) {
			if((value + x) > MAX) {
				value -= x;
				asc   = false;
			} else
				value += x;
		} else {
			if((value - x) < MIN) {
				value += x;
				asc   = true;
			} else
				value -= x;
		}
		temp[i] = value;
		
		//RESULT
		result[i] = 0.0;
	}
}

void usage(int argc, char **argv) {
	std::cerr << "Usage:" << argv[0] << "<order> <iterations>\n";
	std::cerr << "\t<order>      - order for the grid - size= <order>X<order> (positive integer)\n";
	std::cerr << "\t<iterations> - number of iterations\n";
    std::cerr << "\t<pyramid_height> - pyramid heigh\n";
	std::abort();
}

int main(int argc, char** argv) {
    printf("WG size of kernel = %d X %d\n", BLOCK_SIZE, BLOCK_SIZE);

    int order_size, iterations;
    int pyramid_height = 1; // number of iterations
	FLOAT *temp, *power, *result;
	
	/* check validity of inputs*/
	if(argc != 4)
		usage(argc, argv);
	if( (order_size = atoi(argv[1])) <= 0 ||
		(iterations = atoi(argv[2])) <= 0 ||
        (pyramid_height = atoi(argv[3])) <= 0)  
		usage(argc, argv);
    
    /* --------------- pyramid parameters --------------- */
    # define EXPAND_RATE 2// add one iteration will extend the pyramid base by 2 per each borderline
    int borderCols = (pyramid_height)*EXPAND_RATE/2;
    int borderRows = (pyramid_height)*EXPAND_RATE/2;
    int smallBlockCol = BLOCK_SIZE-(pyramid_height)*EXPAND_RATE;
    int smallBlockRow = BLOCK_SIZE-(pyramid_height)*EXPAND_RATE;
    int blockCols = order_size/smallBlockCol+((order_size%smallBlockCol==0)?0:1);
    int blockRows = order_size/smallBlockRow+((order_size%smallBlockRow==0)?0:1);

    printf("pyramidHeight: %d\ngridSize: [%d, %d]\nborder:[%d, %d]\nblockGrid:[%d, %d]\ntargetBlock:[%d, %d]\n",\
	pyramid_height, order_size, order_size, borderCols, borderRows, blockCols, blockRows, smallBlockCol, smallBlockRow);

    /* allocate memory for the temperature and power arrays	*/
	temp   = new FLOAT[order_size * order_size];
	power  = new FLOAT[order_size * order_size];
	result = new FLOAT[order_size * order_size];
	if(!temp || !power || ! result)
		fatal("unable to allocate memory");
	
	/* generate input */
	generate_input(temp, power, result, order_size);

    float *MatrixTemp[2], *MatrixPower;
    cudaMalloc((void**)&MatrixTemp[0], sizeof(FLOAT) * order_size * order_size);
    cudaMalloc((void**)&MatrixTemp[1], sizeof(FLOAT) * order_size * order_size);
    cudaMalloc((void**)&MatrixPower, sizeof(FLOAT) * order_size * order_size);

    long long start_time = get_time();
    cudaMemcpy(MatrixTemp[0], temp, sizeof(FLOAT) * order_size * order_size, cudaMemcpyHostToDevice);
    cudaMemcpy(MatrixPower, power, sizeof(FLOAT) * order_size * order_size, cudaMemcpyHostToDevice);
    
#ifdef VERBOSE
	printf("Start computing the transient temperature\n");
#endif
    
    int ret = compute_tran_temp(MatrixPower, MatrixTemp, order_size, order_size, \
	 iterations, pyramid_height, blockCols, blockRows, borderCols, borderRows);
	printf("Ending simulation\n");
    cudaMemcpy(result, MatrixTemp[ret], sizeof(FLOAT)* order_size * order_size, cudaMemcpyDeviceToHost);

    long long end_time = get_time();
  
    float total_time = ((float) (end_time - start_time)) / (1000*1000);

    //printf("BENCH=hotspot;backend=%s;size=%d;block_size=%d;iterations=%d;threads=%d;gpus=%d;time=%.8f\n", 
	printf("hotspot;%s;%d;%d;%d;%d;%d;%.8f\n",
		"CUDA", order_size, order_size, iterations, 0, 1, total_time);

    cudaFree(MatrixPower);
    cudaFree(MatrixTemp[0]);
    cudaFree(MatrixTemp[1]);
    
    delete[] temp;
	delete[] power;
	delete[] result;
	
	return 0;
}
