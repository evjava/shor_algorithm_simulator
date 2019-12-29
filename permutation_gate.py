from utils import *
from constant_matrices import *

E2SP = csr_matrix(E2)
perm_cache = {
    1: E2SP,
    2: csr_matrix(SWAP_MATRIX)
}

SW = SWAP_SPARSE


def build_perm_matrix(n):
    tms = lambda i: i
    if n % 2 == 0:
        swnm1 = power_kron(SW, n // 2 - 1)
        swn = krons(SW, swnm1)
        prt = krons(E2SP, swnm1, E2SP)
        dot_parts = dots(swn, prt)
        res = power_dot(dot_parts, n // 2)
    else:
        n1 = (n - 1) // 2
        swn = power_kron(SW, n1)
        prt = krons(E2SP, swn)
        prt2 = krons(swn, E2SP)
        dot_parts = dots(prt, prt2)
        res = dots(power_dot(dot_parts, n1), prt)
    return res


def perm_matrix(n):
    res = perm_cache.get(n)
    if res is None:
        res = perm_cache[n] = build_perm_matrix(n)
    return res


def permutation_matrix(t_a, t_b):
    perm_mat_a = perm_matrix(t_a)
    perm_mat = krons(perm_mat_a, sparse.identity(2 ** t_b))
    return perm_mat
