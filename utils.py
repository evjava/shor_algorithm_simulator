import math
import scipy as sp
import numpy as np

from scipy import sparse
from scipy.sparse import csr_matrix
from scipy.stats import rv_discrete


def log2(x):
    return math.log(x, 2)


norm = sp.linalg.norm

def normalize(state):
    return state / norm(state)


def flatten(nested_list):
    return (k for j in nested_list for k in j)

def lflatten(nested_list):
    return list(flatten(nested_list))

def lmap(func, *iterables):
    return list(map(func, *iterables))


def choose_from_probs(probs):
    return rv_discrete(values=(range(len(probs)), probs)).rvs(size=1)[0]


def total_prob(state):
    return np.abs(state ** 2).sum()


def zuka(x):
    return np.vectorize(lambda i: np.round(i, 2))(x)


def almost_zero(x, n=15):
    if n <= 0:
        raise KeyError('positive n expected..')
    return abs(x) < 10 ** (-n)


def replace_empty(x):
    return x if x else '.'


def is_sparse(mat):
    return not isinstance(mat, np.ndarray)


def to_sparse_if_needed(prefix, array):
    if not is_sparse(array):
        # print(prefix, array.shape)
        array = csr_matrix(array)
    return array


def operation(op_name, op, op_sparse):
    def inner(*elements):
        # todo reduce?
        res = None
        for elem in elements:
            if res is None:
                res = elem
            else:
                if isinstance(res, np.ndarray) and isinstance(elem, np.ndarray):
                    res = op(res, elem)
                else:
                    elem = to_sparse_if_needed(f'here elem({op_name})!!!', elem)
                    res = to_sparse_if_needed(f'here res({op_name})!!!', res)
                    try:
                        res = op_sparse(res, elem)
                    except ValueError:
                        raise KeyError(f'Something goes wrong! shapes: {res.shape}, {elem.shape}')
        return res

    return inner


dots = operation('dot', np.dot, csr_matrix.dot)
krons = operation('kron', np.kron, sparse.kron)


def power_op(op):
    def inner(x, times):
        res = x
        for i in range(times - 1):
            res = op(res, x)
        return res

    return inner


power_dot = power_op(dots)
power_kron = power_op(krons)

import time
def timeit(method):
    def timed(*args, **kwargs):
        t_start = time.time()
        result = method(*args, **kwargs)
        t_end = time.time()
        print('\t' * 10, '{}  {:2.2f} s'.format(method.__name__, (t_end - t_start)))
        return result
    return timed
