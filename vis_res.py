#!/usr/bin/env python

import matplotlib.pyplot as plt
import numpy as np
import pylab as pl
from matplotlib import collections  as mc
import sys

lines=[]
x = []
y = []

with open(sys.argv[1], 'r') as f:
    for line in f:
        data = map(float, line.split())
        #print(data)
        lines.append([(data[0], data[1]), (data[2], data[3])])
        x.append(data[0])
        x.append(data[2])
        y.append(data[1])
        y.append(data[3])
        
lc = mc.LineCollection(lines)
fig, ax = pl.subplots()
ax.add_collection(lc)
ax.autoscale()
ax.margins(0.1)
#pl.show()
pl.scatter(x,y)

intersections_x = []
intersections_y = []
with open(sys.argv[2], 'r') as f:
    for line in f:
        data = map(float, line.split())
        intersections_x.append(data[0])
        intersections_y.append(data[1])

pl.scatter(intersections_x, intersections_y, color='r')
pl.show()
#draw points

