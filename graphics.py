#!/usr/bin/env python
import math
import numpy as np
import matplotlib.pyplot as plt

N = 41
hs = (1 - 0)/(N-1)
args = []
vals = []
for i in range(0,N):
	args.append(0 + i*hs)
	vals.append(4*args[i]*(1-args[i]))
data = np.loadtxt("vals_test.txt", delimiter="\n")

print(vals, data)
plt.plot(args, vals, 'r', args, data, 'b')
plt.title('Red line is accuracy solution (4x(1-x)) and blue line is approximate')
plt.show()


