#!/usr/bin/python3.11

import numpy as np
from random import random, seed

#Create hermitian matrix
def sysmake(size, _min, _max):
	mat = [[None for i in range(0, size)] for j in range(0, size)]
	
	for i in range(0, size):
		mat[i][i] = _min + random()*(_max-_min)
		for j in range(i+1, size):
			re = _min + random()*(_max-_min)
			im = _min + random()*(_max-_min)
			mat[i][j] = re + 1j*im
			mat[j][i] = re - 1j*im

	return mat

#M[n][n] -> M[n^2]
def acc2(n, i, j):
	return i*n + j

SIZE = 6
MIN = -10.
MAX = +10.
seed(9642069)

A = sysmake(SIZE, MIN, MAX)
va, ve = np.linalg.eig(A)

print("i, j, ACC2(i,j), A[i][j]:")
for i in range(0, SIZE):
	for j in range(0, SIZE):
		print(i, j, acc2(SIZE, i, j), A[i][j])

for i in range(0, SIZE):
	for j in range(0, SIZE):
		p = '+' if A[i][j].real > 0 else ''
		n = '+' if A[i][j].imag > 0 else ''
		print("M[" + str(acc2(SIZE, i, j)) + "].re = (fp)" + p + str(A[i][j].real) +"; " +\
		      "M[" + str(acc2(SIZE, i, j)) + "].im = (fp)" + n + str(A[i][j].imag) +";")
		

print("\neigvals:")
for i in range(0, SIZE):
	print(va[i])

print("\neigvecs:")
for i in range(0, SIZE):
	print(i, ve[:,i])

print("\nerrs:")
for i in range(0, SIZE):
	print(np.dot(A, ve[:,i]) - va[i]*ve[:,i])




