import numpy as np
import itertools
import random

def createPopulation(n, Nrot, Nworkers):
	"""
	:param n: number of individuals that comprises the population
				 Nrot: number of rotation shifts
				 Nworkers: number of workers and workplaces
	:return: population as an array of matrixes
	"""
	pop = {}
	for i in range(n):
		pop[n] = dict(ind=createIndividual(Nrot, Nworkers), score=0)

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
	:return: matrix that won the tournamet
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


#main code-------------------------------------------

Nworkers = 12
Nshifts = 4

Nmuscls = 1


pop = createPopulation(10, 4, 12)
#calculate fitness of that population



# cost of work position
cstW = 10*np.random.rand(Nworkers)
