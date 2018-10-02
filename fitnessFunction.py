import numpy as np

"""
Derive a function that can be used to measure the fitness of a sequence of shifts.
The function has to be aware that:

 - each shift has a weight for its duration and position on the sequence
 - each working position has a weight representing the effort
 - a sequence with a high effort followed by an high effort should be penalized

D - values for the duration
P - values for the position
E - values for the effort
TA - values that evaluate the sequence of effort
we - weights of E
wta - weights of TA

The calculation of the weight of sequence of effort is derived from the derivative.

The final formula is the following:

exponential((sum(DP*(we*E)) + wta*TA))/1000)

On this example:

[1, 3, 1, 2] - 38
[7, 1, 6, 2] - 347
[7, 8, 9, 1] - 48533
[7, 6, 1, 2] - 735

"""

def weE(seq, we):
	return seq*we

def wtaTA(seq, wta):
	# dseq = np.diff(seq)/(np.diff(seq))
	# for i in range(np.shape(dseq)[1]):
	# 	l = dseq[:, i]*((seq[:, i+1] + seq[:, i]))
	# 	dseq[:, i] = l
	dseq = np.diff(seq)
	print("matrix of diff with sum----------------------------------------------")
	print(dseq)
	dseq = np.hstack([np.zeros((4, 3)), dseq])[:, -4:]
	print("...with ones----------------------------------------------")
	print(dseq)
	# dseq = np.exp(dseq)
	# print("...exponential----------------------------------------------")
	# print(dseq)

	return wta*dseq

def E_TA(seq_E, seq_TA):
	return seq_E + seq_TA

def fitness2():
	t = (DP*we_e)+wta_ta
	return np.sum(t, axis=1).astype(dtype=np.int)

def fitnessMatrix(e_ta, dp):
	# a = np.round(dp * e_ta / 10, 2).astype(dtype=np.int)
	a = dp*e_ta
	print(a.astype(dtype=np.int))
	return np.sum(a, axis=1).astype(dtype=np.int)



DP = [2, 6, 4.5, 6]

Eff = 10*np.random.rand(0, 10)

We = 100
Wta = 10

seqs = np.array([[1, 3, 1, 2], [7,1, 6, 2], [7, 8, 9, 1], [7, 6, 1, 2]])

we_e = weE(seqs, We)
wta_ta = wtaTA(seqs, Wta)
print("matrix of TA")
print(wta_ta)
e_ta = E_TA(we_e, wta_ta)
print("matrix of ETA")
# print(e_ta)
print(fitness2())
# ftns = fitnessMatrix(e_ta, DP)
# print(ftns)
print(np.exp(fitness2()/1000).astype(dtype=np.int))