#!/usr/bin/env python

import sys
import uteis
import argparse
from datetime import datetime

#Variaveis globais
plot_file = "\t\'"+"../scripts/stats_hotspot_speedup.csv"+"\'"
plot_er_margin = " using 1:(${0}-${1}):(${0}+${1}) with filledcu lc rgb color_em notitle, ''"
plot_linha = " using 1:{0} with linespoints t columnheader {1},\\" 

def hotspot_calc_stats(dados, tipo="none", n=0, nb=0, threads=0, gpus=0, tempo_serial=-1):
	tipo_exec = 0
	size = 1
	nbsize = 2
	nthreads = 4
	ngpus = 5
	res = 6
	linhas = list(filter(lambda x: tipo in x[tipo_exec] and str(n) in x[size] and str(nb) in x[nbsize] and str(threads) == x[nthreads] and str(gpus) == x[ngpus], dados))
	if len(linhas) < 1:
		return 0.0, 0.0, 0.0
	amostra = []
	for a in linhas:
		#amostra.append(float(a[res].split(':').pop(1)))
		if tempo_serial > 0:
			amostra.append(tempo_serial / float(a[res])) #Speedup 
			#if gpus == 1 or threads == 40:
			#	print(tempo_serial / float(a[res]))
		else:
			amostra.append(float(a[res]))
	
	return uteis.calc_stats(amostra)

def hotspot_speedup_cpu(dados):
	print("Gerando hotspot speedup...")
	
	nsizes = [12288]
	nblocks = [12288]
	parallel   = "OPENMP"
	sequential = "OPENMP"
	ncpus = [2, 4, 8, 12, 16, 20, 40]

	saida = open("stats_hotspot_speedup.csv", 'w')
	
	#Header e Plotagem
	saida.write("CPUs")
	plot_index = 2
	plot_index_ls = 1
	for n in nsizes:
		for nb in nblocks:
			saida.write(";OpenMP;sd;error")
			print(plot_file \
				+ plot_er_margin.format(plot_index, plot_index+2)\
				+ plot_linha.format(plot_index, "ls " + str(plot_index_ls)))
			plot_index += 3
			plot_index_ls += 1
	
	#Tempo serial
	tempo_serial, sd, error = hotspot_calc_stats(dados, sequential, 12288, 12288, 1, 0)
	print(tempo_serial, sd, error) 
	
	#Dados do arquivo
	saida.write("\n")
	for t in ncpus:
		saida.write(str(t))
		for n in nsizes:
			for nb in nblocks:
				#print(calc_stats(dados, sequential, n, nb, 1))
				mean, sd, error = hotspot_calc_stats(dados, parallel, n, nb, t, 0, tempo_serial)
				saida.write(";" + str(mean) + ";" + str(sd) + ";" + str(error))
		saida.write("\n")

	saida.close()
	print("Concluido \n")

def hotspot_speedup_cuda(dados):
	print("Gerando hotspot cuda...")
	
	nsizes = [12288]
	nblocks = [12288]
	parallel   = "CUDA"
	sequential = "OPENMP"

	saida = open("stats_hotspot_speedup_cuda.csv", 'w')
	
	#Header e Plotagem
	saida.write("Idx")
	for n in nsizes:
		for nb in nblocks:
			saida.write(";GPU (CUDA);sd;error")
	
	#Tempo serial
	tempo_serial, sd, error = hotspot_calc_stats(dados, sequential, 12288, 12288, 1, 0)
	print(tempo_serial, sd, error) 
	
	#Dados do arquivo
	saida.write("\n")
	saida.write("30")
	for n in nsizes:
		for nb in nblocks:
			#print(calc_stats(dados, sequential, n, nb, 1))
			mean, sd, error = hotspot_calc_stats(dados, parallel, n, nb, 0, 1, tempo_serial)
			saida.write(";" + str(mean) + ";" + str(sd) + ";" + str(error))
	saida.write("\n")

	saida.close()
	print("Concluido \n")


def main():
	dados = [] # array das linhas
	parser = argparse.ArgumentParser(description='Estatistica.')
	parser.add_argument("-e", "--entrada", nargs='+', help="arquivos de entrada")
	args = parser.parse_args()
	uteis.le_csv_desempenho(args.entrada, "hotspot", dados)
	
	#print(dados)
	
	hotspot_speedup_cpu(dados)
	hotspot_speedup_cuda(dados)
	
	
if __name__ == "__main__":
	main()

