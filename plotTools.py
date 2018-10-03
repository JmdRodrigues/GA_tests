import numpy as np
import matplotlib.patches as patches
import matplotlib.collections as coll
# import seaborn
# from matplotlib.mlab import find
# import pandas as pd
# import regex as re
import matplotlib.pyplot as plt
# from mpl_toolkits.axes_grid.axes_grid import AxesGrid
from mpl_toolkits.axes_grid.anchored_artists import AnchoredText
from mpl_toolkits.axes_grid1 import make_axes_locatable
# from circles import circles
# from matplotlib.colors import ListedColormap

def make_rgb_transparent(rgb, bg_rgb, alpha):
	return [alpha * c1 + (1 - alpha) * c2
			for (c1, c2) in zip(rgb, bg_rgb)]

def arrayMultiply(array, c):
	return [element*c for element in array]

def arraySum(a, b):
	return map(sum, zip(a,b))

def intermediate(a, b, ratio):
	aComponent = arrayMultiply(a, ratio)
	bComponent = arrayMultiply(b, 1-ratio)
	sumS = [x + y for x, y in zip(aComponent, bComponent)]
	return (sumS[0]/255,sumS[1]/255, sumS[2]/255)

def gradient(a, b, steps):
	steps = [n/float(steps) for n in range(steps)]
	colormap = []
	for step in steps:
		colormap.append(intermediate(a, b, step))
	return colormap

def findcloser(array, value):
	idx = (np.abs(array-value)).argmin()
	return idx


def createSquareMatrix(ncols, nrows, post_colors, force_colors, stations, ind):

	inbetween = 0.1
	wid = 1
	hei = 1
	xx = np.arange(0, ncols, (wid + inbetween))
	yy = np.arange(0, nrows, (hei + inbetween))


	ax = plt.subplot(111, aspect='equal')
	for xi, x in zip(xx, range(0, ncols)):
		for yi, y in zip(yy, range(0, nrows)):
			sq = patches.Rectangle((xi, yi), wid, hei, fc='lightsteelblue', alpha=0.5)
			cr_post = patches.Circle((xi+wid / 4, yi+hei / 2), wid/6, fc=post_colors[y][x], alpha=0.5)
			cr_force = patches.Circle((xi + 3*wid / 4, yi + hei / 2), wid / 6, fc=force_colors[y][x],alpha=0.7)
			ax.add_patch(sq)
			ax.add_patch(cr_post)
			ax.add_patch(cr_force)
			ax.text(xi+wid/2-wid/10, yi+3*hei/4, stations[ind[y][x]-1])

	plt.axis('off')
	ax.axis([0, ncols+1,0,nrows+1])
	plt.show()





