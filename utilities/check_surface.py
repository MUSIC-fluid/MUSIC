#!/usr/bin/env python

count_err = 0
count_lines = 0
f_regulated = open('surface_regulated.dat', 'w')
with open('surface.dat') as f:
    for line in f:
        count_lines += 1
        a = line.split(' ')
        if len(a) == 28:
            f_regulated.write(line)
        else:
            count_err += 1
print(float(count_err)/float(count_lines))
f_regulated.close()
