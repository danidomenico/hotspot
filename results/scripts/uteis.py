#!/usr/bin/env python

import sys
import math
import numpy
from scipy.stats import t
import csv

def le_csv_desempenho(arquivos, prefixo, array_dados):
	if len(arquivos) < 1:
		return
	for a in arquivos:
		with open(a) as f:
			reader = csv.reader(f, delimiter=';', quoting=csv.QUOTE_NONE)
			filtro = list(filter(lambda x: len(x) > 0 and x[0] == prefixo, reader))
			
			for linha in filtro:
				array_dados.append(linha[1:])

def calc_stats(amostra):
	# confidence interval of 95%
	tdist = t.ppf(0.95, len(amostra)-1)
	mean = numpy.mean(amostra)
	std = numpy.std(amostra)
	error = tdist*(std/math.sqrt(len(amostra)))
	return mean, std, error

def calc_speedup_ideal(threads):
	return float(threads)

