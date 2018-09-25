import numpy as np
import itertools
import random

def skiprow(mat1, mat2, Nworkers):
	child1 = mat1
	child2 = mat2

	divider = Nworkers//4

	#mating process--------------------------
	#child1--------------
	child1[:divider, :] = mat2[:divider, :]
	child1[divider:3*divider, :] = mat2[divider:3*divider, :]
	child1[3*divider:, :] = mat2[3*divider:, :]
	#child2--------------
	child2[:divider, :] = mat1[:divider, :]
	child2[divider:3 * divider, :] = mat1[divider:3 * divider, :]
	child2[3 * divider:, :] = mat1[3 * divider:, :]

	return child1, child2

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

def FitnessFunction():


def crossover(ind1, ind2, type):
	#mating step
	if type is "skip_row":
		"""
			With this each row is sequentially shifted between the two individuals
		"""



	elif type is "point_row":


def check_uniqueness(mat, r):
	uniq_r = np.unique(r)
	# all values
	all_v = np.array([i for i in range(1, len(mat) + 1)])
	# possible values
	poss_r = [i for i in all_v if i not in uniq_r]

	return poss_r

def errorCorrection(mat):

	for i in range(len(mat)):
		r = mat[i, :]

		# unique values
		uniq_r = np.unique(r)
		#all values
		all_v = np.array([i for i in range(1, len(mat)+1)])
		# possible values
		poss_r = [i for i in all_v if i not in uniq_r]

		#if the length of unique rotations is not the same as the number of rotations
		if(len() != len(r)):
			#count elements and save index:
			cnt = [i for i in range(len(r)) if sum(r==r[i])>1]

			if(len(cnt) == 2):
				#select one of them randomly...
				ind = cnt[random.randint(0, 1)]
				#...and change it:
				mat[i, :][ind] = poss_r[random.randint(0, len(poss_r)-1)] #minus 1 because the randint function is inclusive


def TournamentSelection():

def createPopulation(n, Nrot, Nworkers):
	"""

	:param n: number of individuals that comprises the population
				 Nrot: number of rotation shifts
				 Nworkers: number of workers and workplaces
	:return: population as an array of matrixes
	"""
	for i in range(n):
		createIndividual(Nrot, Nworkers)


# def inidividual():

Nworkers = 12
Nshifts = 4
Nmuscls = 1



# cost of work position
cstW = 10*np.random.rand(Nworkers)



# print(cstW)



