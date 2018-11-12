import numpy as np
import matplotlib.patches as patches
from decimal import Decimal
import matplotlib.collections as coll
# import seaborn
# from matplotlib.mlab import find
# import pandas as pd
# import regex as re
import matplotlib.pyplot as plt
import seaborn as sns
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

def add_squared_matrix(ax, xx, yy, ncols, nrows, wid, hei, all_colors, ind, stations, median_list, workers, wp_dict, Circles):
	# initial colors for diversity
	bls = sns.color_palette('Blues')
	green = (60/255, 179/255, 113/255)
	gold = (1,215/255,0)
	orange = 'orange'
	red = (178/255, 34/255, 34/255)
	for xi, x in zip(xx, range(0, ncols)):
		ax.text(xi + wid/3, nrows + hei + 0.5, "Rot " + str(x))
		for yi, y in zip(reversed(yy), range(0, nrows)):
			sq_exp = patches.Rectangle((xi, yi + hei / 2), wid, hei / 2,
									   fc=all_colors["Score_total_estacao"][ind[y][x] - 1], alpha=1)
			ax.add_patch(sq_exp)
			ax.text(xi + wid / 2 - len(str(wp_dict[stations[ind[y][x] - 1]]['OS'])) / 10, yi + 3 * hei / 5, str(wp_dict[stations[ind[y][x] - 1]]['OS']))
			tk = wp_dict[stations[ind[y][x] - 1]]['%B']+wp_dict[stations[ind[y][x] - 1]]['%TB']
			sh = wp_dict[stations[ind[y][x] - 1]]['%SH']+wp_dict[stations[ind[y][x] - 1]]['%HL']
			el = wp_dict[stations[ind[y][x] - 1]]['%R6']

			scores_div = [tk, sh, el]
			br_div = ['tk','sh', 'el']

			med = median_list[workers[y]]
			for xi2 in [0, 1, 2]:
				div_x = xi + xi2 * wid / 3
				if (scores_div[xi2] >= med[br_div[xi2]]):
					color = red
				else:
					color = green

				if(xi2<2):
					sq_div = patches.Rectangle((div_x, yi), (wid / 3) - wid / 100, hei / 2,
										   fc=color, alpha=0.3)
					ax.add_patch(sq_div)
					ax.text(div_x + wid / 6 - len(str(round(Decimal(scores_div[xi2]), 2))) / 10, yi + hei / 10, str(round(Decimal(scores_div[xi2]), 2)))

			sq_div = patches.Rectangle((xi + 2 * wid / 3, yi), (wid / 3), hei / 2,
									   fc=color, alpha=0.3)
			ax.text(xi + 2 * (wid / 3) + (wid / 6) - len(str(el)) / 10, yi + hei / 10, str(el))
			ax.add_patch(sq_div)
			if (Circles == True):
				for index, key in enumerate(all_colors.keys()):
					cr_post = patches.Circle((xi + 0.5 + 1.3 * index / 5, yi + hei / 2), 1 / 7,
											 fc=all_colors[key][ind[y][x] - 1], alpha=0.3)
					ax.add_patch(cr_post)

			ax.text(-1.2, yi + hei / 2, "Wkr" + str(y))
# ax.text(-1.2, yi + hei / 2, "Wkr" + str(y))
# ax.text(wid * ncols + 0.2, yi + hei / 2, str(exp_score_per_seq[y]))
# ax.text(wid * ncols + 1.2, yi + hei / 2, str(div_score_per_seq[y]))

def createSquareMatrix2(ncols, nrows, all_colors, stations, ind, exp_score_per_seq, div_score_per_seq, med_list, workers, wp_dict, Circles=False):
	inbetween = 0.1
	wid = 5
	hei = 1
	xx = np.arange(0, wid * ncols, (wid + inbetween))
	yy = np.arange(0, hei * (nrows + 1), (hei + inbetween))
	xx_exp = xx[-1] + inbetween + wid
	xx_exp2 = xx_exp + inbetween +wid/2

	fig = plt.figure(figsize=(100, 100))
	ax = plt.subplot(111, aspect='equal')
	#create matrix of squares
	add_squared_matrix(ax, xx, yy, ncols, nrows, wid, hei, all_colors, ind, stations, med_list, workers, wp_dict, Circles)
	for index, yi in enumerate(reversed(yy)):
		sq_scores_exp = patches.Rectangle((xx_exp, yi+hei/2), wid/2, hei/2, facecolor='grey', alpha=0.7)
		sq_scores_exp2 = patches.Rectangle((xx_exp2, yi), wid/2, hei, facecolor='grey', alpha=0.7)

		ax.add_patch(sq_scores_exp)
		ax.add_patch(sq_scores_exp2)

		ax.text(xx_exp + (wid / 3) - len(str(round(Decimal(exp_score_per_seq[index]),2))) / 10, yi + 2*hei/3, str(round(Decimal(exp_score_per_seq[index]),2)))
		ax.text(xx_exp + wid / 3 - len("Score")/10, nrows + hei + 0.5, "Score")
		worker_final = exp_score_per_seq[index]/(1+div_score_per_seq[index])
		ax.text(xx_exp2 + (wid / 4) - len(str(round(Decimal(worker_final),2))) / 10, yi + hei/2, str(round(Decimal(worker_final),2)))
		ax.text(xx_exp2 + wid / 4 - len("Score C")/10, nrows + hei + 0.5, "Score C")

		sq_scores_div = patches.Rectangle((xx_exp, yi), wid/2, hei/2, facecolor='grey', alpha=0.3)
		ax.add_patch(sq_scores_div)
		print(div_score_per_seq)
		ax.text(xx_exp + (wid / 3) - len(str(round(Decimal(div_score_per_seq[index]),2))) / 10, yi+0.1, str(round(Decimal(div_score_per_seq[index]),2)))
	exp_final_sc = np.sum(exp_score_per_seq)*(1+abs(np.std(exp_score_per_seq)))
	div_final_sc = np.sum(div_score_per_seq)/(1+abs(np.std(div_score_per_seq)))

	ax.text(wid * ncols, -1.0, "Exposure Score: " + str(round(Decimal(exp_final_sc),2)) + "   var: " + str(round(Decimal(abs(np.std(exp_score_per_seq))),2)))
	ax.text(wid * ncols, -1.5, "Diversity Score: " + str(round(Decimal(div_final_sc),2)) + "   var: " + str(round(Decimal(abs(np.std(div_score_per_seq))),2)))
	ax.text(wid * ncols, -0.5, "Final Score: " + str(exp_final_sc/div_final_sc))
	plt.axis('off')
	ax.axis([0, wid*(ncols+2), 0, nrows+2])


def createSquareMatrix(ncols, nrows, all_colors, stations, ind, score_per_seq, Circles=False):
	inbetween = 0.1
	wid = 5
	hei = 1
	xx = np.arange(0, wid*ncols, (wid + inbetween))
	yy = np.arange(0, hei*nrows + 1, (hei + inbetween))

	ax = plt.subplot(111, aspect='equal')
	for xi, x in zip(xx, range(0, ncols)):
		ax.text(xi+wid/4, nrows+hei+1, "Rot "+str(x))
		for yi, y in zip(reversed(yy), range(0, nrows)):
			sq = patches.Rectangle((xi, yi), wid, hei, fc=all_colors["Score_total_estacao"][ind[y][x]-1], alpha=1)
			ax.add_patch(sq)
			if(Circles == True):
				for index, key in enumerate(all_colors.keys()):
					cr_post = patches.Circle((xi + 0.5 + 1.3*index/5, yi+hei / 2), 1/7, fc=all_colors[key][ind[y][x]-1], alpha=0.5)
					ax.add_patch(cr_post)
				# cr_post2 = patches.Circle((xi + 2 / 5, yi + hei / 2), 1 / 10, fc=all_colors["%_Bent"][ind[y][x]-1], alpha=0.5)
			# cr_post3 = patches.Circle((xi + 3 / 5, yi + hei / 2), 1 / 10, fc=all_colors["%_twisted_lateral_bent"][ind[y][x]-1], alpha=0.5)
			# cr_post4 = patches.Circle((xi + 4 / 5, yi + hei / 2), 1 / 10, fc=all_colors["%Shoulder height"][ind[y][x]-1], alpha=0.5)
			# cr_post5 = patches.Circle((xi + 5 / 5, yi + hei / 2), 1 / 10, fc=all_colors["%_head_level"][ind[y][x]-1], alpha=0.5)
			# cr_post6 = patches.Circle((xi + 6 / 5, yi + hei / 2), 1 / 10, fc=all_colors["%_range60"][ind[y][x]-1], alpha=0.5)
			# cr_post7 = patches.Circle((xi + 7 / 5, yi + hei / 2), 1 / 10, fc=all_colors["%_range80"][ind[y][x]-1], alpha=0.5)
			# cr_post8 = patches.Circle((xi + 8 / 5, yi + hei / 2), 1 / 10, fc=all_colors["%_range100"][ind[y][x]-1], alpha=0.5)
			# cr_post9 = patches.Circle((xi + 9 / 5, yi + hei / 2), 1 / 10, fc=all_colors["Force_Score"][ind[y][x]-1], alpha=0.5)
			# cr_post10 = patches.Circle((xi + 10 / 5, yi + hei / 2), 1 / 10, fc=all_colors["Force_light"][ind[y][x]-1], alpha=0.5)
			# cr_post11 = patches.Circle((xi + 11 / 5, yi + hei / 2), 1 / 10, fc=all_colors["Force_midddle"][ind[y][x]-1], alpha=0.5)
			# cr_post12 = patches.Circle((xi + 12 / 5, yi + hei / 2), 1 / 10, fc=all_colors["Force_hard"][ind[y][x]-1], alpha=0.5)
			# cr_post13= patches.Circle((xi + 13 / 5, yi + hei / 2), 1 / 10, fc=all_colors["Force_stressed"][ind[y][x]-1], alpha=0.5)
			# cr_post14 = patches.Circle((xi + 14 / 5, yi + hei / 2), 1 / 10, fc=all_colors["Score_vibrations"][ind[y][x]-1], alpha=0.5)


			# ax.add_patch(cr_post2)
			# ax.add_patch(cr_post6)
			# ax.add_patch(cr_post9)
			# ax.add_patch(cr_post14)
			ax.text(xi+wid/2-wid/6, yi+2*hei/5, stations[ind[y][x]-1])
			ax.text(-1.8,yi+hei/2, "Wkr" + str(y))
			ax.text(wid*ncols+1.2, yi+hei/2, str(score_per_seq[y]))
		ax.text(wid*ncols, -0.8, "Final Score:" + str(np.sum(score_per_seq)))

	plt.axis('off')
	ax.axis([0, wid*ncols+1,0,nrows+1])
	# plt.show()

def createWorkplaceMatrix(risk_factors):

	inbetween = 0.1
	wid = 2
	hei = 1
	ncols = len(list(risk_factors.keys())[3:])
	nrows = len(risk_factors["URQ"])
	stations = risk_factors["Estação"]
	risk_factors_tags = ["P_S", "%B", "%SB", "%TB", "%SH", "%HL", "%R6", "%R8", "%RT", "FS", "FL", "FM", "FH", "FS", "SV", "LS", "OS"]
	xx = np.arange(0, 2*ncols, (wid + inbetween))
	yy = np.arange(0, nrows, (hei + inbetween))
	ax = plt.subplot(111, aspect='equal')

	colors_factors = {}

	tag = 0
	for key, xi in zip(list(risk_factors.keys())[4:], xx):
		col = convert2color(key, risk_factors[key])
		colors_factors[key] = col
		for y, yi in zip(range(0, len(col)), yy):
			sq = patches.Rectangle((xi, yi), wid, hei, fc=col[y], alpha=0.5)
			ax.add_patch(sq)
			# add text of workplace
		ax.text(xi+wid/8, 12 * hei + 1, risk_factors_tags[tag], size=14)
		tag+=1


	#for the overall score:
	col = convert2color("Score_total_estacao", risk_factors["Score_total_estacao"])
	colors_factors["Score_total_estacao"] = col
	for y, yi in zip(range(0, len(col)), yy):
		sq_w = patches.Rectangle((max(xx) + 1, yi), wid, hei, fc="white", alpha=0.5)
		ax.add_patch(sq_w)
		sq = patches.Rectangle((max(xx)+2, yi), wid, hei, fc=col[y], alpha=0.5)
		ax.add_patch(sq)
		ax.text(-1.8 * wid, yi + hei / 4, stations[y])
	ax.text(wid*ncols + 2, 12 * hei + 1, risk_factors_tags[-1], size=14)

	plt.axis('off')
	ax.axis([0, 2*(ncols+3), 0, nrows+2])
	plt.show()

	return colors_factors


def convert2color(key, col):
	green = (60, 179, 113)
	gold = (255,215,0)
	red = (178, 34, 34)

	if(key == "Score_total_estacao"):
		max_s = 80
		min_s = 30

		ll3_5 = np.linspace(30, 50, 50000)
		color1 = gradient(red, gold, 50000)
		color2 = gradient(gold, green, 50000)
		color = color2 + color1
		# print(len(color))
		colors = []
		for s in col:
			if(s < 30):
				colors.append(color2[0])
			elif(s >30 and s<50):
				colors.append(color1[findcloser(ll3_5, s)])
			else:
				colors.append(color1[-1])
	else:
		if(key == "Postura_score"):
			max_s = 50
			min_s = 0
		elif(key == "Force_Score"):
			max_s = 50
			min_s = 0
		elif(key == "Score_vibrations"):
			max_s = 4
			min_s = 0
		elif("%" in key):
			max_s = 100
			min_s = 0
		else:
			min_s = 0
			max_s = 1

		ll = np.linspace(min_s, max_s, 100000)
		color1 = gradient(red, gold, 50000)
		color2 = gradient(gold, green, 50000)
		color = color2+color1
		# print(len(color))

		colors = [color[findcloser(ll, s)] for s in col]

	return colors