#!/usr/bin/env python

import sys
import uteis
import argparse
from datetime import datetime

#Variaveis globais
plot_file = "\t\'"+"../scripts/stats_ideal_speedup.csv"+"\'"
plot_linha = " using 1:{0} with linespoints t columnheader {1},\\" 

def gera_speedup_file(ncpus):
	print("Gerando speedup ideal...")
	
	saida = open("stats_ideal_speedup.csv", 'w')
	
	#Header e Plotagem
	saida.write("CPUs")
	
	plot_index = 2
	plot_index_ls = 13
	
	#Dados cabecalho e plot
	saida.write(";\"" + str("Ideal") + "\"")
	print(plot_file \
		 + plot_linha.format(plot_index, "ls " + str(plot_index_ls)))
	
	#Dados do arquivo
	saida.write("\n")
	for t in ncpus:
		saida.write(str(t))
		#Speedup ideal
		mean = uteis.calc_speedup_ideal(t)
		saida.write(";" + str(mean))
		saida.write("\n")

	saida.close()
	print("Concluido \n")


def main():
	dados = [] # array das linhas
	parser = argparse.ArgumentParser(description='Estatistica.')
	parser.add_argument("-t", "--threads", nargs='+', help="threads para arquivo")
	args = parser.parse_args()
	
	gera_speedup_file(args.threads)
	
	
if __name__ == "__main__":
	main()

