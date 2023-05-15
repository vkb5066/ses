#!/usr/bin/python3.11

from matplotlib import pyplot as plt

c, f = [], []
with open("data.tst", 'r') as infile:
	for lin in infile:
		line = lin.split()
		if(len(line) != 6):
			continue 
		c.append(float(line[2][:-1]))
		f.append(float(line[5]))
	infile.close()

plt.plot(c, f, 'k.')
plt.show()
