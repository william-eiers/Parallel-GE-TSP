import math
import sys
import numpy as np

def dist(x1,y1,x2,y2):
	return math.sqrt(float((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1)))

num = 0
cities = []
with open('wi29.tsp') as f:
	for line in f:
		data = line.split()
		if len(data) == 1:
			num = int(data[0])
		else:
			cities.append((float(data[1]),float(data[2])))

distances = [[0 for x in range(num)] for x in range(num)]
target = open('distances3','w')
target.write(str(num))
for x in range(0,num):
	target.write("\n")
	for y in range(0,num):
		if x == y:
			distances[x][y] = 0
		else:
			distances[x][y] = dist(cities[x][0], cities[x][1], cities[y][0],cities[y][1])
		target.write(str(distances[x][y]))
		target.write(" ")
target.close()
