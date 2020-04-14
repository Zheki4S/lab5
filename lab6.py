import numpy as np
import math
def fx0(fx, i0, x, h):
	n = len(fx[1])
	matrix = [[0] * n for i in range(n)];
	y = [0 for i in range(n)]
	y[0] = 3 * (fx[1][1] - fx[1][0]);
	for i in range(1, n - 1):
		y[i] = 3 * (fx[1][i + 1] - fx[1][i - 1])
	y[n - 1] = 3 * (fx[1][n - 1] - fx[1][n - 2])
	for i in range(n):
		for j in range(n):
			if(i == j):
				matrix[i][j] = 4
			if(i == j + 1 or i == j - 1):
				matrix[i][j] = 1
	matrix[0][0] = 2
	matrix[n - 1][n - 1] = 2				
	k = np.linalg.solve(matrix, y)
	n = 4
	a0 = [[1, 0, 0, 0], [0, 1, 0, 0], [1, 1, 1, 1], [0, 1, 2, 3]];
	y0 = [fx[1][i0], k[i0], fx[1][i0 + 1], k[i0 + 1]]
	a = np.linalg.solve(a0, y0)
	f = 0
	for i in range(4):
		f += a[i] * (x - i0 * h) ** i
	return f;
def trap(fx, h0, k, N):
	trap = 0
	h = k * h0
	for i in range(0, N - 1, k):
		trap += (fx[1][i] + fx[1][i + k])
	return h * trap / 2
def simp(fx, h0, k, N):
	simp = 0
	h = k * h0
	for i in range(0, N - 1, 2 * k):
		simp += fx[1][i] + 4 * fx[1][i + k] + fx[1][i + 2 * k]
	return h * simp / 3
def m3d8(fx, h0, k, N):
	m3d8 = 0
	h = k * h0
	for i in range(0, N - 1, 3 * k):
		m3d8 += fx[1][i] + 3 * (fx[1][i + k] + fx[1][i + 2 * k]) + fx[1][i + 3 * k]	
	return h * m3d8 * 3 / 8	
def runge(foo, fx, h0, N):
	return	math.fabs(foo(fx, h0, 1, N) - foo(fx, h0, 2, N))
		
N = 13
h0 = 0.125
fx = [[0, 0.125, 0.25, 0.375, 0.5, 0.625, 0.75, 0.875, 1], [0, 0.02147, 0.29305, 0.494105, 0.541341, 0.516855, 0.468617, 0.416531, 0.367879]]
nanoh = 0.00001
h1 = 1 / (N - 1)
fxN = [[], []]
k = 0
for i in range(len(fx[1]) - 1):
	for j in range(1, int(h0 / nanoh) + 1):
		x = i * h0 + j * nanoh
		if(math.fabs(x - k * h1) < 2 * nanoh):
			fxN[0].append(x)
			fxN[1].append(fx0(fx, i, x, h0))
			k += 1
print("trap: " + str(trap(fxN, h1, 1, N)))
print("runge trap: " + str(runge(trap, fxN, h1, N) / 3))
print("simp: " + str(simp(fxN, h1, 1, N)))
print("runge simp: " + str(runge(simp, fxN, h1, N) / 15))	
print("3 / 8: " + str(m3d8(fxN, h1, 1, N)))
print("runge 3 / 8: " + str(runge(m3d8, fxN, h1, N) / 31))	


		