#!/usr/bin/env python

from numpy import *
from os import path

raw_file_path = "../initial/u_field_1.dat"

raw_file = open(raw_file_path, 'r')
header = raw_file.readline()
nx = int(header.split()[6])
ny = int(header.split()[8])
dx = float(header.split()[12])
dy = float(header.split()[14])
print("raw grid size: %d, raw grid spacing: %.6f" % (nx, dx))

grid_out_nx = 301
grid_shift = int((grid_out_nx - 1 - nx)/2)
e_out = zeros([grid_out_nx, grid_out_nx]) + 1e-12
utau_out = zeros([grid_out_nx, grid_out_nx]) + 1.
ux_out = zeros([grid_out_nx, grid_out_nx])
uy_out = zeros([grid_out_nx, grid_out_nx])
x_min = y_min = 0.
for i in range(nx):
    for j in range(ny):
        data = [float(x) for x in raw_file.readline().split()]
        if i == 0 and j == 0:
            x_min = data[1]
            y_min = data[2]
        e_out[grid_shift+i, grid_shift+j] = data[3]
        utau_out[grid_shift+i, grid_shift+j] = data[4]
        ux_out[grid_shift+i, grid_shift+j] = data[5]
        uy_out[grid_shift+i, grid_shift+j] = data[6]
    raw_file.readline()
raw_file.close()

# output file
output_file_path = "../initial/u_field_1_enlarged.dat"
output_file = open(output_file_path, "w")
header = ("#0.0 0.0 0.0 n_eta= 4 nx= %d ny= %d deta= 5.0 dx= %.6f dy= %.6f\n"
          % (grid_out_nx, grid_out_nx, dx, dy))
output_file.write(header)
for i in range(grid_out_nx):
    x_local = x_min + (i - grid_shift)*dx
    for j in range(grid_out_nx):
        y_local = y_min + (j - grid_shift)*dy
        output_file.write(
            "0.0  %.8e  %.8e  %.8e  %.8e  %.8e  %.8e  %.8e  %.8e  %.8e  %.8e\n"
            % (x_local, y_local, e_out[i, j], utau_out[i, j], ux_out[i, j],
               uy_out[i, j], 0.0, 0.0, 0.0, 0.0))
output_file.close()
