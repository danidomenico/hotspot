# C compiler
CC = g++
#ICC = icc
#CC_FLAGS = -fopenmp -O3
CC_FLAGS = -g -O3 -std=c++11 -fopenmp
#OFFLOAD_CC_FLAGS = -offload-option,mic,compiler,"-no-opt-prefetch"

#all: hotspot hotspot_offload
all: hotspot


hotspot: hotspot_openmp.cpp Makefile 
	$(CC) $(CC_FLAGS) hotspot_openmp.cpp -o hotspot.openmp 

#hotspot_offload: hotspot_openmp.cpp Makefile
#	$(ICC) $(CC_FLAGS) $(OFFLOAD_CC_FLAGS) -DOMP_OFFLOAD hotspot_openmp.cpp -o hotspot_offload

clean:
	rm -f hotspot.openmp
