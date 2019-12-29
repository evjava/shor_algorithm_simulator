import numpy as np

from scipy.sparse import csr_matrix

HADAMARD_MATRIX = np.array([[1, 1], [1, -1]]) / np.sqrt(2)
E2 = np.identity(2)
HM_SPARSE = csr_matrix(HADAMARD_MATRIX)
SWAP_MATRIX = np.array([[1, 0, 0, 0], [0, 0, 1, 0], [0, 1, 0, 0], [0, 0, 0, 1]])
SWAP_SPARSE = csr_matrix(SWAP_MATRIX)
CNOT = np.array([[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 0, 1], [0, 0, 1, 0]])
