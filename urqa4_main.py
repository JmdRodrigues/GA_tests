import pandas as pd
import numpy as np
from GA_tests.plotTools import gradient, findcloser, createSquareMatrix
import matplotlib.pyplot as plt
import random

"""

1 - load exposure data from apergo and team leader opinion
2 - calculate a fitness value for an example of sequences

"""
def fitnessSeq(wp_seq, wp_dict):
	keys = ['p_s', 'f_s']
	seq_temp = []
	for key in keys:
		seq_k = np.array([wp_dict[i][key] for i in wp_seq])
		seq_temp.append(fitnessCalc(seq_k))
	return np.sum(seq_temp)

def fitnessCalc(seq):
	DP = [1, 4, 2, 3]
	We = 100

	Wta = 10

	we_e = We * seq
	wta_ta = variability(seq, Wta)

	t = np.sum((DP * we_e) + wta_ta)

	return t

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


def crossover(ind1, ind2, nworkers, type):
	#mating step
	if type is "skip_row":
		"""
			With this each row is sequentially shifted between the two individuals
		"""
		skiprow(ind1, ind2, nworkers)

	# elif type is "point_row":

def TournamentSelection(pop, k):
	"""
	:param k: tournament size
	:param pop: entire population
	:return: matrix that won the tournament
	"""

	if k > len(pop):
		print("the tournament pool size should be lower than the population size")
	else:
		tournamentSel = random.sample(range(0, len(pop)), k)

		tournamentPop = pop[tournamentSel]

		#Compete the elements of the population
		# winner = the one with the higher score

def skiprow(mat1, mat2, Nworkers):

	child1 = mat1
	child2 = mat2
	divider = Nworkers // 4
	# mating process--------------------------
	# child1--------------
	child1[:divider, :] = mat2[:divider, :]
	child1[divider:3 * divider, :] = mat2[divider:3 * divider, :]
	child1[3 * divider:, :] = mat2[3 * divider:, :]
	# child2--------------
	child2[:divider, :] = mat1[:divider, :]
	child2[divider:3 * divider, :] = mat1[divider:3 * divider, :]
	child2[3 * divider:, :] = mat1[3 * divider:, :]

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
	wp_dict = {}
	stations = apergo["Estação"].to_dict()
	print(stations)
	posture_score = apergo["Postura_score"].to_dict()
	force_score = apergo["Force_Score"].to_dict()

	for i in range(0, len(stations)):
		wp_dict[stations[i]] = {'p_s': posture_score[i],
								'f_s': force_score[i]
								}
	return wp_dict, stations, posture_score, force_score

def createScoreMatrixANDColor(score, ind):
	"""

	:param score: score dictionary
	:param ind: individual, represented as a matrix with sequences in each row
	:return: matrix with all the scores
	"""
	scr_mat = []
	clr_mat = []

	green = (60, 179, 113)
	red = (178, 34, 34)

	ll = np.linspace(min(score.values()), max(score.values()), 100000)
	color = gradient(red, green, 100000)

	print(score)
	for seq in ind:
		seq_i =  [score[ws-1] for ws in seq]
		scr_mat.append(seq_i)
		#color correspondance
		colors = [color[findcloser(ll, s)] for s in seq_i]
		clr_mat.append(colors)

	return scr_mat, clr_mat


wp_dict, stations, post_score, force_score = load_scores()

ind = createIndividual(4, 12)
post_mat, post_color = createScoreMatrixANDColor(post_score, ind)
force_mat, force_color = createScoreMatrixANDColor(force_score, ind)

print(post_mat)
print(post_color[0][0])

createSquareMatrix(4, 12, post_color, force_color, stations, ind)


# sts = [stations[i] for i in seq]
# sts2 = [stations[i] for i in seq2]
# sts3 = [stations[i] for i in seq3]
#
# print(fitnessSeq(sts, wp_dict))
# print(fitnessSeq(sts2, wp_dict))
# print(fitnessSeq(sts3, wp_dict))

