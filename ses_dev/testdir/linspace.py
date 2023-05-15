#!/usr/bin/python3.11

from numpy import linspace

#--- 1d -------------------------------------------
BEG, END, N = 1.2, 2.6, 11

def li1(b, e, n):
	ret = [None for i in range(0, n)]

	dx = (e - b)/(n - 1)
	for i in range(0, n):
		ret[i] = b + dx*i

	return ret
"""	
nu = linspace(BEG, END, N)
mi = li1(BEG, END, N)
print("from", BEG, "to", END, "in", N)
print("numpy:\n")
print(nu)
print("mine:\n")
print(mi)
"""


#--- 2d -----------------------------------------
BEG, END, N = [5.6, 1.4], [1.2, 2.8], 6

def li2(b, e, n):
	ret = [[None, None] for i in range(0, n)]

	dx = (e[0] - b[0])/(n - 1)
	dy = (e[1] - b[1])/(n - 1)
	for i in range(0, n):
		ret[i][0] = b[0] + dx*i
		ret[i][1] = b[1] + dy*i

	return ret
	
nu = linspace(BEG, END, N)
mi = li2(BEG, END, N)
print("from", BEG, "to", END, "in", N)
print("numpy:\n")
print(nu)
print("mine:\n")
print(mi)


