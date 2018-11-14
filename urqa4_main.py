import pandas as pd
import numpy as np
from GA_tests.plotTools import gradient, findcloser, createSquareMatrix, createWorkplaceMatrix, createSquareMatrix2, linear_plot_quantile
import matplotlib.pyplot as plt
import random
import time

"""

1 - load exposure data from apergo and team leader opinion
2 - calculate a fitness value for an example of sequences

"""

"""
Added new fitness calculation. The purpose is to convert the current model to a multi-objective function
which comprises:
- maximizing diversity
- minimizing exposure

Diversity is the extent at which exposure entities differ. In this case, we want to evaluate 
diversity as the change of exposure levels in body regions.
For this, we should be aware that:

- it is independent of the exposure level
- it should only depend on the successive changes in diversity

Questions:
- is going from 0 to 47 a good level of diversity?
- is going from 0 to 10 worse or better than the previous one?
- is there any difference that accounts for diversity?

I believe that in this case, we should look at transitions of levels to check diversity.
We can start by defining the mean +/- the standard deviation and get 4 levels, which give us
the linear distribution of exposure values for that body region on all workstations.
After that, any level transition is given a score of 1, else, it is 0 and there is no diversity.

Next, we have to check the homogenety in diversity for the all body...this should be calculateed
based on weights...because shoulder exposure might be worse...at least, we may check if the ap ergo
gives weights to body regions...if so, we could use these weights.
The final equation for homogenety can be:
-->  (Ws*shoulder_div + Wt*trunk_div + We*elbow_div)/std (std is dividing, because higher values of std, mean more dispersion,
which should be penalized)

By minimizing exposure, we can reduce the chance of haveing sequences with high exposures. In this case, the exposure will be calculated
with the final score given by the ap ergo function. 
"""

def objective_exposure(exposure_mat, stations, diversity_d, diversity_mstd):
	scores_exp = []
	scores_div = []
	scores_div_per_musc = {'tk':[], 'sh':[], 'el':[]}
	for worker in exposure_mat:
		seq = [stations[w - 1] for w in worker]
		div_seq_scores, div_per_musc = objective_diversity(seq, diversity_d, diversity_mstd)
		for key in div_per_musc.keys():
			scores_div_per_musc[key].append(div_per_musc[key])
		scores_div.append(np.sum(div_seq_scores))
		scores_exp.append(exposure_calc(seq))

	f_sdiv = np.sum(scores_div)/(1+np.std(scores_div))
	f_sexp = np.sum(scores_exp)*(1+np.std(scores_exp))
	print(f_sdiv)
	#get standard deviation for homogeneity of scores
	return scores_exp, f_sexp, scores_div, scores_div_per_musc, f_sdiv

def exposure_calc(exposure_seq):
	duration = [0.23, 0.38, 0.19, 0.20]
	seq_k = np.array([wp_dict[i]["OS"] for i in exposure_seq])

	return np.sum(np.multiply(duration, seq_k))

def objective_exposureV3(exposure_mat, stations, diversity_d, diversity_quant, workers):
	scores_exp = []
	scores_div = []
	scores_div_per_musc = {'tk': [], 'sh': [], 'el': []}
	quant_sec_per_workr = []
	for index, worker in enumerate(exposure_mat):
		seq = [stations[w - 1] for w in worker]
		div_final_score, div_per_br, quant_sec = objective_diversityV4(workers[index], seq, diversity_d, diversity_quant)
		for key in div_per_br.keys():
			scores_div_per_musc[key].append(div_per_br[key])
		scores_div.append(div_final_score)
		scores_exp.append(exposure_calc(seq))
		quant_sec_per_workr.append(quant_sec)

	f_sdiv = np.sum(scores_div) / (1 + abs(np.std(scores_div)))
	f_sexp = np.sum(scores_exp) * (1 + abs(np.std(scores_exp)))

	# get standard deviation for homogeneity of scores
	return scores_exp, f_sexp, scores_div, scores_div_per_musc, f_sdiv, quant_sec_per_workr

def objective_exposureV2(exposure_mat, stations, diversity_d, diversity_median):
	scores_exp = []
	scores_div = []
	scores_div_per_musc = {'tk':[], 'sh':[], 'el':[]}
	for worker in exposure_mat:
		seq = [stations[w - 1] for w in worker]
		div_seq_scores, div_per_musc = objective_diversityV3(seq, diversity_d, diversity_median)
		for key in div_per_musc.keys():
			scores_div_per_musc[key].append(div_per_musc[key])
		scores_div.append(np.sum(div_seq_scores))
		scores_exp.append(exposure_calc(seq))

	f_sdiv = np.sum(scores_div)/(1+np.std(scores_div))
	f_sexp = np.sum(scores_exp)*(1+np.std(scores_exp))
	print(f_sdiv)
	#get standard deviation for homogeneity of scores
	return scores_exp, f_sexp, scores_div, scores_div_per_musc, f_sdiv


# def objective_diversityV2(seq, diversity_d, diversity_median):
# 	for ws in seq:
# 		for br in diversity_d.keys():
# 			if(dd[ws][br] >= mstd[br][0]):
# 				ss[br].append(4)
# 			elif(dd[ws][br] < mstd[br][0] and dd[ws][br] > mstd[br][1]):
# 				ss[br].append(3)
# 			elif(dd[ws][br] <mstd[br][1] and dd[ws][br] > mstd[br][2]):
# 				ss[br].append(2)
# 			elif (dd[ws][br] < mstd[br][2]):
# 				ss[br].append(1)
# 	distance_med = divers

def objective_diversity(exposure_seq, diversity_d, diversity_mstd):
	"""
	divide in 4 levels of diversity (mean +/- std)
	:return:
	"""
	ss = get_diversity_levels(exposure_seq, diversity_d, diversity_mstd)
	dss = []
	for i in ss.keys():
		dss.append(len(np.where(np.diff(ss[i])!=0)[0]))

	return dss/(1+np.std(dss)), ss

def get_diversity_levels(seq, dd, mstd):
	ss = {'tk':[], 'sh':[], 'el':[]}
	for ws in seq:
		for br in mstd.keys():
			if(dd[ws][br] >= mstd[br][0]):
				ss[br].append(4)
			elif(dd[ws][br] < mstd[br][0] and dd[ws][br] > mstd[br][1]):
				ss[br].append(3)
			elif(dd[ws][br] <mstd[br][1] and dd[ws][br] > mstd[br][2]):
				ss[br].append(2)
			elif (dd[ws][br] < mstd[br][2]):
				ss[br].append(1)

	return ss

def objective_diversityV4(worker, seq, dd, quantile):
	# score for each body region
	f_ss = {'tk': [], 'sh': [], 'el': []}
	quant_sec_per_br = {'tk': [], 'sh': [], 'el': []}
	for br in ["tk", "sh", "el"]:
		quant_seq = []
		for ws in seq:
			if (dd[ws][br] >= quantile[worker][br]['75']):
				quant_seq.append(4)
			elif (dd[ws][br] >= quantile[worker][br]['50']) and (dd[ws][br] < quantile[worker][br]['75']):
				quant_seq.append(3)
			elif (dd[ws][br] >= quantile[worker][br]['25']) and (dd[ws][br] < quantile[worker][br]['50']):
				quant_seq.append(2)
			else:
				quant_seq.append(0)

		# diversity for each body region (not counting with the amplitude of the jump
		quant_sec_per_br[br] = quant_seq
		f_ss[br] = sum([1 if(df>0) else 0 for df in abs(np.diff(quant_seq))])

	final_score = (f_ss["tk"] + f_ss["sh"] + f_ss["el"]) / (1 + np.std(f_ss["tk"] + f_ss["sh"] + f_ss["el"]) ** 2)

	return final_score, f_ss, quant_sec_per_br


def objective_diversityV3(worker, seq, dd, medd):
	#score for each body region
	f_ss = {'tk':[], 'sh':[], 'el':[]}

	for br in ["tk", "sh", "el"]:
		med_seq = []
		for ws in seq:
			if(dd[ws][br] >= medd[worker][br]):
				med_seq.append(1)
			else:
				med_seq.append(0)

		# diversity for each body region
		f_ss[br] = sum(abs(np.diff(med_seq)))

	final_score = (f_ss["tk"] + f_ss["sh"] + f_ss["el"])/(1+np.std(f_ss["tk"] + f_ss["sh"] + f_ss["el"])**2)

	return final_score, f_ss


def fitness(ind, stations, wp_dict):
	val = []
	for worker in ind:
		seq = [stations[w-1] for w in worker]
		score = fitnessSeq(seq, wp_dict)
		val.append(score)
	# print(val)
	ind_score = sum(val)
	return val, ind_score

def fitnessSeq(wp_seq, wp_dict):
	keys = ["P_S", "%B", "%SB", "%TB", "%SH", "%HL", "%R6", "%R8", "%RT", "FS", "FL", "FM", "FH", "FS", "SV", "LS", "OS"]
	weight = [10, 0, 0, 0, 0, 0, 0, 0, 0, 10, 0, 0, 0, 0, 10, 0, 0]
	seq_temp = []
	for key in keys:
		seq_k = np.array([wp_dict[i][key] for i in wp_seq])
		seq_temp.append(fitnessCalc(seq_k))
	return np.sum(np.multiply(weight, seq_temp))

def fitnessCalc(seq):
	cumulation = [1, 3, 2, 4]
	duration = [0.23, 0.38, 0.19, 0.20]

	DP = np.multiply(cumulation, duration)
	# DP = duration
	We = 100

	Wta = 10

	we_e = We * seq
	wta_ta = variability(seq, Wta)

	t = np.sum((DP*we_e) + wta_ta)

	return np.exp(t/100)

def variability(seq, w):
	dseq = np.diff(seq)
	dseq = np.insert(dseq, 0, 0)

	return w * dseq

def createPopulation(n, Nrot, Nworkers):
	"""
	:param n: number of individuals that comprises the population
				 Nrot: number of rotation shifts
				 Nworkers: number of workers and workplaces
	:return: population as an array of matrixes
	"""
	pop = {}
	for i in range(n):
		pop[i] = dict(ind=createIndividual(Nrot, Nworkers), score=0)

	return pop

def createIndividual(Nshifts, Nworkers):
	"""Make one attempt to generate a filled m**2 x m**2 Sudoku board,
	returning the board if successful, or None if not.
	"""

	board = np.array([[None for _ in range(Nshifts)] for _ in range(Nworkers)])
	for i in range(Nshifts):
		numbers = random.sample(range(1, Nworkers + 1), Nworkers)

		if (i == 0):
			board[:, 0] = numbers
		else:
			notValid = True
			while (notValid):
				notValid = False
				for ii in range(Nworkers):
					existing_values = board[ii, :]
					if (numbers[ii] in existing_values):
						notValid = True
						numbers = random.sample(range(1, Nworkers + 1), Nworkers)

			board[:, i] = numbers
	return board

# def FitnessFunction():


def crossover(ind1, ind2, nworkers, nshifts):
	types = ["skip_row", "skipcol"]
	# tp = types[random.randint(0, len(types)-1)]
	tp = random.randint(1, 10)
	print(tp)
	#mating step
	if tp <= 3:
		"""
			With this each row is sequentially shifted between the two individuals
		"""
		mat1, mat2 = skiprow(ind1, ind2, nworkers)
		mat1 = GeneControl(mat1)
		mat2 = GeneControl(mat2)
		d = 0

	elif tp > 3:
		mat1, mat2 = skipcol(ind1, ind2, nshifts)
		mat1 = GeneControl(mat1)
		mat2 = GeneControl(mat2)

		d = 1

	return mat1, mat2, d


def TournamentSelection(pop, k, n):
	"""
	:param k: tournament size
	:param n: number of winners
	:param pop: entire population
	:return: matrix that won the tournament
	"""

	if k > len(pop):
		print("the tournament pool size should be lower than the population size")
	else:
		winners = [0, 0]
		while(np.multiply(winners[0], winners[1]) != 0):
			tournamentSel = random.sample(range(0, len(pop)), k)
			winners = [0, 0]
			previous_score = 0
			for round in range(0, n):
				min_score = 100000000
				for candidate in tournamentSel:
					ind_score = pop[candidate]['score']
					if(np.logical_and(ind_score < min_score, ind_score > previous_score)):
						min_score = ind_score
						winners[round] = candidate
						previous_score = min_score

				tournamentSel.remove(winners[round])

		return winners

		#Compete the elements of the population
		# winner = the one with the higher score

def skipcol(prnt1, prnt2, Nshifts):
	"""
	not ready for other nshifts than 4
	:param prnt1:
	:param prnt2:
	:param Nshifts:
	:return:
	"""

	prnt1 = np.array(prnt1["ind"])
	prnt2 = np.array(prnt2["ind"])

	child1 = prnt1
	child2 = prnt2

	rnd1_1 = random.randint(1, Nshifts)
	rnd1_2 = random.randint(1, Nshifts)
	rnd2_1 = random.randint(1, Nshifts)
	rnd2_2 = random.randint(1, Nshifts)

	lst = range(0, Nshifts)
	node1 = rnd1_1-1
	node2 = rnd2_1-1
	lst_n1 = np.r_[lst[node1:], lst[:node1]].astype(int)
	lst_n2 = np.r_[lst[node2:], lst[:node2]].astype(int)

	if(rnd1_2%2 == 0):
		#pair means that it is alternate
		if(rnd2_2%2 == 0):
			#pair means that it is alternate
			for i in range(0, Nshifts, 2):
				print(lst_n1[i])
				print(lst_n2[i])
				child1[:, lst_n1[i]] = prnt2[:, lst_n2[i]]
				child2[:, lst_n1[i]] = prnt1[:, lst_n2[i]]
		else:
			#even means that it is sequential
			for i, ii in zip(range(0, Nshifts, 2), range(0, Nshifts//2)):
				child1[:, lst_n1[i]] = prnt2[:, lst_n2[ii]]
				child2[:, lst_n1[i]] = prnt1[:, lst_n2[ii]]

	else:
		# pair means that it is alternate
		if (rnd2_2 % 2 == 0):
			# pair means that it is alternate
			for i, ii in zip(range(0, Nshifts, 2), range(0, Nshifts//2)):
				child1[:, lst_n1[ii]] = prnt2[:, lst_n2[i]]
				child2[:, lst_n1[ii]] = prnt1[:, lst_n2[i]]
		else:
			# even means that it is sequential
			for i in range(0, Nshifts // 2):
				child1[:, lst_n1[i]] = prnt2[:, lst_n2[i]]
				child2[:, lst_n1[i]] = prnt1[:, lst_n2[i]]

	return child1, child2

def skiprow(mat1, mat2, Nworkers):

	child1 = mat1["ind"]
	child2 = mat2["ind"]
	divider = Nworkers // 4
	# mating process--------------------------
	# child1--------------
	child1[:divider][:] = child2[:divider][:]
	child1[divider:3 * divider][:] = child2[divider:3 * divider][:]
	child1[3 * divider:][:] = child2[3 * divider:][:]
	# child2--------------
	child2[:divider][:] = mat1['ind'][:divider][:]
	child2[divider:3 * divider][:] = mat1['ind'][divider:3 * divider][:]
	child2[3 * divider:][:] = mat1['ind'][3 * divider:][:]

	return child1, child2

# def check_uniqueness(mat, r):
# 	uniq_r = np.unique(r)
# 	# all values
# 	all_v = np.array([i for i in range(1, len(mat) + 1)])
# 	# possible values
# 	poss_r = [i for i in all_v if i not in uniq_r]
#
# 	return poss_r

def GeneControl(mat):

	a = mat
	# go through each row to find repeated values
	for i in range(np.shape(a)[0]):
		r = a[i, :]
		# unique elements of the row
		uniq_r = np.unique(r)
		# all possible values
		all_v = np.array([ii for ii in range(1, len(a) + 1)])
		# possible values for that row
		poss_r = [ii for ii in all_v if ii not in uniq_r]

		# if size of unique vector is not the same size as the row, it is because
		# there is a repeated value
		if (len(uniq_r) != len(r)):
			# index of repeated values
			ind = [ii for ii in range(len(r)) if sum(r == r[ii]) > 1]

			# entire column with first repeated value
			col = a[:, ind[0]]
			# repeated value to be swapped
			val = r[ind[0]]

			# for ii in each possible value
			for ii in poss_r:
				# row and col where there is a possible swapping action
				row_ind = np.where(col == ii)[0]
				col_ind = ind[0]

				# row where the ii value is (this row is to find out if the repeated value can be put in here)
				rr = a[np.where(col == ii)[0], :]

				# can the repeated value be put in this row. If so...
				if (val not in rr):
					a[i, col_ind] = ii
					a[row_ind, col_ind] = val
					break

	return mat

def load_scores(): #load apergo data - Fatores de risco AP ergo
	apergo = pd.read_excel("URQA4/Fatores_risco_APErgo_URQA4.xlsx")

	stations = apergo["Estação"].to_dict()

	posture_score = apergo["Postura_score"].to_dict()
	force_score = apergo["Force_Score"].to_dict()
	vibration_score = apergo["Score_vibrations"].to_dict()
	key_tags = ["OS", "P_S", "%B", "%SB", "%TB", "%SH", "%HL", "%R6", "%R8", "%RT", "FS", "FL", "FM", "FH", "FS", "SV", "LS"]

	wp_dict = {station:{tag:0 for tag in key_tags} for station in stations.values()}

	tag_dict = {}

	#normalize vectors
	for tag, key in zip(key_tags, apergo.keys()[3:]):
		mn = apergo[key].mean()
		std = apergo[key].std()
		score_list = np.array(list(apergo[key]))
		# if(std == 0):
		# 	score_list = np.array(list(apergo[key]))
		# else:
		# 	score_list = (np.array(list(apergo[key])) - mn)/std
		tag_dict[tag] = score_list


	for i in range(0, len(stations)):
		for tag in key_tags:
			wp_dict[stations[i]][tag] = tag_dict[tag][i]


	return wp_dict, stations, posture_score, force_score, vibration_score, apergo

def createScoreMatrixANDColor(score, ind):
	"""
	:param score: score dictionary
	:param ind: individual, represented as a matrix with sequences in each row
	:return: matrix with all the scores
	"""
	scr_mat = []
	clr_mat = []

	green = (60, 179, 113)
	gold = (255, 215, 0)
	red = (178, 34, 34)

	ll = np.linspace(min(score.values()), max(score.values()), 100000)
	color1 = gradient(red, gold, 50000)
	color2 = gradient(gold, green, 50000)
	color = color2+color1

	for seq in ind:
		seq_i =  [score[ws-1] for ws in seq]
		scr_mat.append(seq_i)
		#color correspondance
		colors = [color[findcloser(ll, s)] for s in seq_i]
		clr_mat.append(colors)

	return scr_mat, clr_mat

def fitnesspop(pop, stations, wp_dict):
	for index in range(0, len(pop)):
		# seq_scores, ind_score = fitness(pop[index]['ind'], stations, wp_dict)
		# print(seq_scores)
		exp_score_per_seq, exp_ind_score, div_score_per_seq, div_scores_per_sec_musc, div_ind_score = objective_exposure(
			pop[index]['ind'], stations, diversity_dict, diversity_mstd)
		final_score = (np.sum(exp_score_per_seq) * np.std(exp_score_per_seq)) / (
					np.sum(div_score_per_seq) / np.std(div_score_per_seq))
		pop[index]['score'] = final_score

	return pop

def newpop(mat1, mat1_score, mat2, mat2_score, poplt):
	max_array = []
	keys = list(poplt.keys())
	n = 0
	while(n < 2):
		max_score=0
		for index, key in enumerate(keys):
			ind_score = poplt[key]['score']
			if(ind_score > max_score):
				max_score = ind_score
				maxi = index

		max_array.append(keys[maxi])
		keys.remove(keys[maxi])
		n+=1

	if mat2_score < poplt[max_array[1]]['score']:
		poplt[max_array[1]].update({'ind':mat2, 'score':mat2_score})
		poplt[max_array[0]].update({'ind': mat1, 'score': mat1_score})
	else:
		if mat1_score < poplt[max_array[1]]['score']:
			poplt[max_array[1]].update({'ind': mat1, 'score': mat1_score})

			if(mat2_score < poplt[max_array[0]]['score']):
				poplt[max_array[0]].update({'ind': mat2, 'score': mat2_score})

	return poplt

def meanScores(pop):
	sum_s = 0
	for key in pop.keys():
		sum_s += pop[key]['score']

	return sum_s/len(pop)

def get_best_pop(poplt):
	keys = list(pop.keys())
	max_score = 0
	for index, key in enumerate(keys):
		ind_score = poplt[key]['score']
		if (ind_score > max_score):
			max_score = ind_score
			maxi = index
	return keys[maxi]

def get_best_pop2(poplt):
	keys = list(pop.keys())
	max_score = 100000
	for index, key in enumerate(keys):
		ind_score = poplt[key]['score']
		if (ind_score < max_score):
			max_score = ind_score
			maxi = index
	return keys[maxi]

def load_TL_example():
	# example matrix from Team Leader

	code_dict = dict(A=1, B=3, C=4, D=2, E=5, F=6, H=7, I=8, J=9, K=10, L=11,
					 M=12)

	tl_matrix = [['A', 'K', 'B', 'I'],
				 ['E', 'D', 'J', 'A'],
				 ['C', 'J', 'E', 'H'],
				 ['B', 'M', 'L', 'E'],
				 ['K', 'B', 'M', 'F'],
				 ['J', 'C', 'H', 'M'],
				 ['D', 'I', 'K', 'B'],
				 ['M', 'L', 'F', 'D'],
				 ['L', 'H', 'C', 'J'],
				 ['H', 'F', 'A', 'K'],
				 ['I', 'E', 'D', 'C'],
				 ['F', 'A', 'I', 'L'], ]
	tl_matrix2 = [
		["K",	"D", "C", "H"],
				  ["C",	"H", "I", "E"],
				  ["I",	"E", "D", "L"],
				  ["E",	"J", "A", "F"],
				  ["D", "M", "E", "K"],
				  ["M", "I", "F", "C"],
				  ["L", "K", "M", "M"],
				  ["A", "F", "B", "A"],
				  ["H", "A", "J", "B"],
				  ["J", "C", "L", "D"],
				  ["B", "L", "K", "I"],
				  ["F", "B", "H", "J"]
	              ]

	tl_matrix3 = [["M", "H", "A", "M"],
	              ["A", "L", "J", "E"],
	              ["I", "A", "H", "C"],
	              ["F", "B", "I", "J"],
	              ["B", "F", "C", "H"],
	              ["J", "M", "L", "D"],
	              ["K", "I", "M", "B"],
	              ["C", "D", "E", "K"],
	              ["E", "J", "D", "A"],
	              ["L", "E", "K", "I"],
	              ["D", "K", "B", "F"],
	              ["H", "C", "F", "L"]]

	tl_matrix4 = [["M", "E", "D", "J"],
	              ["D", "I", "E", "M"],
	              ["C", "H", "B", "F"],
	              ["A", "K", "I", "L"],
	              ["H", "C", "J", "K"],
	              ["E", "D", "L", "A"],
	              ["I", "J", "K", "E"],
	              ["B", "M", "C", "D"],
	              ["J", "B", "M", "H"],
	              ["L", "A", "F", "C"],
	              ["F", "L", "H", "B"],
	              ["K", "F", "A", "I"]]

	tl_matrix5 = [
		["L", "B", "K", "D"],
		["H", "L", "J", "E"],
		["M", "D", "D", "K"],
		["E",  "M", "B", "F"],
		["C", "H", "M", "A"],
		["A", "I", "E", "L"],
		["K", "J", "C", "I"],
		["B", "F", "H", "C"],
		["J", "A", "F", "B"],
		["I", "K", "A", "M"],
		["D", "C", "L", "H"],
		["F", "E", "I", "J"]]

	tl_matrx6 = [
		["A", "F", "C", "K"],
		 ["M", "A", "F", "D"],
		["K", "E", "D", "J"],
		["E", "K", "I", "H"],
		["B", "L", "M", "A"],
		["H", "I", "L", "B"],
		["L", "H", "A", "M"],
		["I", "M", "J", "E"],
		["C", "D", "B", "I"],
		["J", "B", "H", "C"],
		["D", "C", "K",  "L"],
		["F", "J", "E", "F"]
	]

	tl_dict_pop = {"1":tl_matrix, '2': tl_matrix2, '3':tl_matrix3, '4':tl_matrix4, '5':tl_matrix5, '6':tl_matrx6}

	for tl in tl_dict_pop.keys():
		ind_tl = [[code_dict[tl_dict_pop[tl][j][i]] for i in range(0, 4)] for j in range(0, 12)]
		tl_dict_pop[tl] = ind_tl

	return tl_dict_pop


def validateMatrix(mat):
	#just checks i matrix is valid (has different line values and different column values)

	#check lines:
	for i in mat:
		if(len(np.unique(i)) == len(i)):
			print("lines are valid")
		else:
			return "not valid"

	#check columns:
	for j in range(0, np.shape(mat)[1]):
		if(len(mat[:, j]) == len(np.unique(mat[:,j]))):
			print("columns are valid")
		else:
			return "not valid"

def list_scores_diversity(score_dict):
	diversity_dict = {}
	diversity_mean_std = {}
	tk = []
	sh = []
	el = []
	for workstation in score_dict.keys():
		diversity_dict[workstation] = {'tk':score_dict[workstation]['%B']+score_dict[workstation]['%TB'],
									   'sh':score_dict[workstation]['%SH']+score_dict[workstation]['%HL'],
									   'el':score_dict[workstation]['%R6']}

		tk.append(score_dict[workstation]['%B']+score_dict[workstation]['%TB'])
		sh.append(score_dict[workstation]['%SH']+score_dict[workstation]['%HL'])
		el.append(score_dict[workstation]['%R6'])

	tk_mean = np.mean(tk)
	tk_std = np.std(tk)
	sh_mean = np.mean(sh)
	sh_std = np.std(sh)
	el_mean = np.mean(sh)
	el_std = np.std(sh)

	diversity_mean_std['tk'] = [tk_mean + tk_std, tk_mean, tk_mean - tk_std]
	diversity_mean_std['sh'] = [sh_mean + sh_std, sh_mean, sh_mean - sh_std]
	diversity_mean_std['el'] = [el_mean + el_std, el_mean, el_mean - el_std]

	return diversity_dict, diversity_mean_std

def diversity_list_2(score_dict):
	diversity_dict = {}
	diversity_mean_std = {}
	tk = []
	sh = []
	el = []
	for workstation in score_dict.keys():
		diversity_dict[workstation] = {'tk': score_dict[workstation]['%B'] + score_dict[workstation]['%TB'],
									   'sh': score_dict[workstation]['%SH'] + score_dict[workstation]['%HL'],
									   'el': score_dict[workstation]['%R6'] + score_dict[workstation]['%R8']
											 + score_dict[workstation]['%RT']}


		tk.append(score_dict[workstation]['%B'] + score_dict[workstation]['%TB'])
		sh.append(score_dict[workstation]['%SH'] + score_dict[workstation]['%HL'])
		el.append(score_dict[workstation]['%R6'] + score_dict[workstation]['%R8']
											 + score_dict[workstation]['%RT'])

	tk_med = np.median(tk)
	sh_med = np.median(sh)
	el_med = np.median(el)

	ax1 = plt.subplot(3,1,1)
	ax1.plot(tk, 'o')
	ax1.set_title("Trunk: "+ str(tk_med))
	ax1.axhline(tk_med, xmin=0, xmax=len(tk))
	ax2 = plt.subplot(3,1,2)
	ax2.plot(sh, 'o')
	ax2.set_title("Shoulder: "+ str(sh_med))
	ax2.axhline(sh_med, xmin=0, xmax=len(sh))
	ax3 = plt.subplot(3,1,3)
	ax3.plot(el, 'o')
	ax3.set_title("Elbow: "+ str(el_med))
	ax3.axhline(el_med, xmin=0, xmax=len(el))
	plt.show()

def load_median_list(score_dict, versatility):
	diversity_dict = {}
	diversity_median = {}

	for workstation in score_dict.keys():
		diversity_dict[workstation] = {'tk': score_dict[workstation]['%B'] + score_dict[workstation]['%TB'],
									   'sh': score_dict[workstation]['%SH'] + score_dict[workstation]['%HL'],
									   'el': score_dict[workstation]['%R6'] + score_dict[workstation]['%R8']
											 + score_dict[workstation]['%RT']}
	for worker in versatility.index:
		med_wk = {br: np.median([diversity_dict[ws][br] for ws in score_dict.keys()]) for br in ['tk', 'sh', 'el']}
		diversity_median[worker] = med_wk

	return diversity_median

def load_percentile_list(score_dict, versatility):
	diversity_dict = {}
	diversity_percentile = {}

	for workstation in score_dict.keys():
		diversity_dict[workstation] = {'tk': score_dict[workstation]['%B'] + score_dict[workstation]['%TB'],
									   'sh': score_dict[workstation]['%SH'] + score_dict[workstation]['%HL'],
									   'el': score_dict[workstation]['%R6'] + score_dict[workstation]['%R8']
											 + score_dict[workstation]['%RT']}
	for worker in versatility.index:
		percentile_wk = {br:{'25': np.percentile([diversity_dict[ws][br] for ws in score_dict.keys()], 25),
						 '50': np.percentile([diversity_dict[ws][br] for ws in score_dict.keys()], 50),
						 '75': np.percentile([diversity_dict[ws][br] for ws in score_dict.keys()], 75)}  for br in ['tk', 'sh', 'el']}

		diversity_percentile[worker] = percentile_wk

	return diversity_percentile

def load_versatility(file):
	xl = pd.read_excel(file, "Versatilidade")
	print(xl)
	return xl.iloc[:12, :12], xl.index[:12]


#.....................................................................
#a fazer para mostrar 3ª:
"""
1 - fazer load de varias matrizes planeadas pelos TL e verificar o resultado atribuido
2 - testar o algoritmo e verificar se gera melhores matrizes
3 - testar o algoritmo com varios tipos de:
				a - crossover
				b - tamanho de população
				c - probabilidades de crossover
				d - mutacoes
				e - de que forma podemos encaixar ai os sintomas
"""

"""Tarefa 1-------------------------------------------------------------------------------------------------------

	Load das matrizes dos team leaders e dar um score
"""
#load scores and versatility matrixes
wp_dict, stations, post_score, force_score, vibration_score, risk_factors = load_scores()
diversity_dict, diversity_mstd = list_scores_diversity(wp_dict)
versatility, workers = load_versatility(r"C:\Users\rjoao\Documents\PhD\Algoritmo Plano Rotacional\GA_tests\URQA4\Matriz versatilidade_A4B.xlsx")
#create list of median list for each worker and body region
median_list = load_median_list(wp_dict, versatility)
#load_percentile list
percentile_list = load_percentile_list(wp_dict, versatility)
print(percentile_list)
#Matrix that shows the color risk of each workplace
all_colors = createWorkplaceMatrix(risk_factors)

#1 - load de matrizes de TL
tl_dct = load_TL_example()

#team leader matrixes representation
for tl in tl_dct.keys():
	exp_score_per_seq, exp_ind_score, div_score_per_seq, div_scores_per_sec_musc, div_ind_score, quant_sec = objective_exposureV3(
		tl_dct[tl], stations, diversity_dict, percentile_list, workers)
	createSquareMatrix2(4, 12, all_colors, stations, tl_dct[tl], exp_score_per_seq, div_score_per_seq, div_scores_per_sec_musc, percentile_list, workers, wp_dict)
	linear_plot_quantile(quant_sec, workers)
	plt.show()
#
# #create a chromossome - rotation matrix
# ind = createIndividual(4, 12)
# #find if it is valid or not
# validateMatrix(ind)

# #calculate exposure and diversity vaues for the chromossome
# exp_score_per_seq, exp_ind_score, div_score_per_seq, div_scores_per_sec_musc, div_ind_score = objective_exposureV3(ind, stations, diversity_dict, median_list, workers)
# #show the results in the matrix
# createSquareMatrix2(4, 12, all_colors, stations, ind, exp_score_per_seq, div_score_per_seq, diversity_mstd, wp_dict)
# plt.show()




#.....................................................................

#genetic algorithm test
pop = createPopulation(50, 4, 12)

for individual in range(0, len(pop)):
	ind = pop[individual]['ind']
	exp_score_per_seq, exp_ind_score, div_score_per_seq, div_scores_per_sec_musc, div_ind_score = objective_exposureV3(
		pop[individual]['ind'], stations, diversity_dict, median_list, workers)
	# createSquareMatrix2(4, 12, all_colors, stations, ind, exp_score_per_seq, div_score_per_seq, diversity_mstd, wp_dict)
	# plt.show()

pop_scores = []
#score each pop
min_score = 100000000
min_index = 0
for index in range(0, len(pop)):
	# seq_scores, ind_score = fitness(pop[index]['ind'], stations, wp_dict)
	exp_score_per_seq, exp_ind_score, div_score_per_seq, div_scores_per_sec_musc, div_ind_score = objective_exposureV3(
		pop[index]['ind'], stations, diversity_dict, median_list, workers)
	# print(seq_scores)
	final_score = (np.sum(exp_score_per_seq)*(1+np.std(exp_score_per_seq)))/(np.sum(div_score_per_seq)/(1+np.std(div_score_per_seq)))
	pop[index]['score'] = final_score
	if(final_score < min_score):
		min_score = exp_ind_score
		min_index = index

# np.savez('10000pops.npz', pop=[pop])
# print("Saved!!!")
# container = np.load("10000pops.npz")
# pop2 = [container[key] for key in container]
# print(pop2)


# post_mat, post_color = createScoreMatrixANDColor(post_score, ind)
# force_mat, force_color = createScoreMatrixANDColor(force_score, ind)

all_colors = createWorkplaceMatrix(risk_factors)

# #create square matrix of the best scored individual
# createSquareMatrix(4, 12, all_colors, stations, pop[min_index]['ind'], score_per_seq)
#
# #team_leader matrix
# createSquareMatrix(4, 12, all_colors, stations, load_TL_example(), score_per_seq_tl)

#create square matrix of the best scored individual
# createSquareMatrix(4, 12, all_colors, stations, pop[min_index]['ind'], score_per_seq)

#team_leader matrix
# createSquareMatrix(4, 12, all_colors, stations, load_TL_example(), score_per_seq_tl)

# plt.show()

worst_pop = get_best_pop(pop)
ind_worse = pop[worst_pop]['ind']

#try the genetic algorithm
n = 1
iterations=1
m_score = [meanScores(pop)]
previous_score = meanScores(pop)
x_data = [0]
x_data2 = [0]
m_score2 = [0]
fig = plt.figure()
ax = fig.add_subplot(211)
ax2 = fig.add_subplot(212)
Ln, = ax.plot(x_data, m_score)
Ln2, = ax2.plot(x_data2, m_score2)
# plt.ion()
nbr_repeatedValues = 0


best_pop1 = get_best_pop2(pop)
best_ind1 = pop[best_pop1]['ind']


while(n<1000):

	#fitness population
	pop = fitnesspop(pop, stations, wp_dict)
	#tournament selection
	winners = TournamentSelection(pop, 5, 2)
	#crossover
	mat1, mat2, d = crossover(pop[winners[0]], pop[winners[1]], 12, 4)
	#fitness of new elements
	exp_score_per_seq, exp_ind_score, div_score_per_seq, div_scores_per_sec_musc, div_ind_score = objective_exposureV3(
		mat1, stations, diversity_dict, median_list, workers)
	# print(seq_scores)
	final_score_mat1 = (np.sum(exp_score_per_seq) * (1 + np.std(exp_score_per_seq))) / (
				np.sum(div_score_per_seq) / (1 + np.std(div_score_per_seq)))
	exp_score_per_seq, exp_ind_score, div_score_per_seq, div_scores_per_sec_musc, div_ind_score = objective_exposureV3(
		mat2, stations, diversity_dict, median_list, workers)
	# print(seq_scores)
	final_score_mat2 = (np.sum(exp_score_per_seq) * (1 + np.std(exp_score_per_seq))) / (
			np.sum(div_score_per_seq) / (1 + np.std(div_score_per_seq)))
	# val_n, mat1_score = fitness(mat1, stations, wp_dict)
	# val_n, mat2_score = fitness(mat2, stations, wp_dict)
	#check if scores are good
	pop = newpop(mat1, final_score_mat1, mat2, final_score_mat2, pop)
	#get mean scores
	m_score.append(meanScores(pop))
	x_data.append(iterations)
	Ln.set_xdata(x_data)
	Ln.set_ydata(m_score)
	Ln2.set_xdata(x_data)
	Ln2.set_ydata(d)
	#increase loopping round
	n+=1
	iterations+=1
	ax.set_xlim(0, iterations+10)
	ax.set_ylim(min(m_score)-10, max(m_score)+10)
	ax2.set_xlim(0, iterations + 10)
	plt.pause(.01)
	plt.xlabel("Nbr of Cycles")
	plt.ylabel("Mean Score of Population")

	# if(m_score[-2]==m_score[-1]):
	# 	nbr_repeatedValues+=1
	# if(nbr_repeatedValues > 20):
	# 	n = 1
	# 	pop = createPopulation(50, 4, 12)
	# 	current_score = meanScores(pop)
	# 	nn = 0
	# 	while(current_score>previous_score or nn<50):
	# 		pop = createPopulation(50, 4, 12)
	# 		current_score = meanScores(pop)
	# 		nn+=1
	# 	previous_score = current_score
	# 	nbr_repeatedValues=0

	if(m_score[-2]==m_score[-1]):
		nbr_repeatedValues+=1
	else:
		nbr_repeatedValues = 0
	if(nbr_repeatedValues > 40):
		n = 1
		pop = createPopulation(50, 4, 12)
		current_score = meanScores(pop)
		nn = 0
		while(current_score>previous_score or nn<10):
			pop = createPopulation(50, 4, 12)
			current_score = meanScores(pop)
			nn+=1
		previous_score = current_score
		nbr_repeatedValues=0

plt.show(block=True)

best_pop2 = get_best_pop2(pop)
best_ind2 = pop[best_pop2]['ind']

print(best_ind1)
print(best_ind2)
score_per_seq_worse, ind_score = fitness(best_ind1, stations, wp_dict)
score_per_seq_best, ind_score2 = fitness(best_ind2, stations, wp_dict)
fig = plt.figure(figsize=(200, 200))
createSquareMatrix(4, 12, all_colors, stations, best_ind1, score_per_seq_worse, subplt=211)
createSquareMatrix(4, 12, all_colors, stations, best_ind2, score_per_seq_best, subplt=212)
plt.suptitle("Initial matrix ("+str(ind_score)+") vs Ending Matrix ("+str(ind_score2)+")")

plt.show(block=True)
