import pandas as pd
import numpy as np
from GA_tests.plotTools import gradient, findcloser, createSquareMatrix, createWorkplaceMatrix
import matplotlib.pyplot as plt
import random
import time

"""

1 - load exposure data from apergo and team leader opinion
2 - calculate a fitness value for an example of sequences

"""
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


def crossover(ind1, ind2, nworkers, type="skip_row"):
	#mating step
	if type is "skip_row":
		"""
			With this each row is sequentially shifted between the two individuals
		"""
		mat1, mat2 = skiprow(ind1, ind2, nworkers)
		mat1 = errorCorrection(mat1)
		mat2 = errorCorrection(mat2)

		return mat1, mat2

	# elif type is "point_row":

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

def errorCorrection(mat):

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
		if(std == 0):
			score_list = np.array(list(apergo[key]))
		else:
			score_list = (np.array(list(apergo[key])) - mn)/std
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
		seq_scores, ind_score = fitness(pop[index]['ind'], stations, wp_dict)
		# print(seq_scores)
		pop[index]['score'] = ind_score


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

	ind_tl = [[code_dict[tl_matrix[j][i]] for i in range(0, 4)] for j in range(0, 12)]

	return ind_tl

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

#.....................................................................
wp_dict, stations, post_score, force_score, vibration_score, risk_factors = load_scores()

# print(wp_dict)

ind = createIndividual(4, 12)
validateMatrix(ind)
score_per_seq, ind_score = fitness(ind, stations, wp_dict)
score_per_seq_tl, ind_score_tl = fitness(load_TL_example(), stations, wp_dict)

#.....................................................................

#genetic algorithm test
pop = createPopulation(10000, 4, 12)

pop_scores = []
#score each pop
min_score = 100000000
min_index = 0
for index in range(0, len(pop)):
	seq_scores, ind_score = fitness(pop[index]['ind'], stations, wp_dict)
	# print(seq_scores)
	pop[index]['score'] = ind_score
	if(ind_score < min_score):
		min_score = ind_score
		min_index = index


np.savez('10000pops.npz', pop=[pop])
print("Saved!!!")
container = np.load("10000pops.npz")
pop2 = [container[key] for key in container]
print(pop2)


# post_mat, post_color = createScoreMatrixANDColor(post_score, ind)
# force_mat, force_color = createScoreMatrixANDColor(force_score, ind)

all_colors = createWorkplaceMatrix(risk_factors)


# #create square matrix of the best scored individual
# createSquareMatrix(4, 12, all_colors, stations, pop[min_index]['ind'], score_per_seq)
#
# #team_leader matrix
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
fig = plt.figure()
ax = fig.add_subplot(111)
Ln, = ax.plot(x_data, m_score)
# plt.ion()
nbr_repeatedValues = 0

best_pop1 = get_best_pop2(pop)
best_ind1 = pop[best_pop1]['ind']

while(n<35):
	#fitness population
	pop = fitnesspop(pop, stations, wp_dict)
	#tournament selection
	winners = TournamentSelection(pop, 5, 2)
	#crossover
	mat1, mat2 = crossover(pop[winners[0]], pop[winners[1]], 12)
	#fitness of new elements
	val_n, mat1_score = fitness(mat1, stations, wp_dict)
	val_n, mat2_score = fitness(mat2, stations, wp_dict)
	#check if scores are good
	pop = newpop(mat1, mat1_score, mat2, mat2_score, pop)
	#get mean scores
	m_score.append(meanScores(pop))
	x_data.append(iterations)
	Ln.set_xdata(x_data)
	Ln.set_ydata(m_score)
	#increase loopping round
	n+=1
	iterations+=1
	ax.set_xlim(0, iterations+10)
	ax.set_ylim(min(m_score)-10, max(m_score)+10)
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

# plt.show()



