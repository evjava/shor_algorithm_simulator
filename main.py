import pdb
from fractions import Fraction
from math import gcd

import fire
from memorize import memorize

from permutation_gate import *
from shor_config import ShorConfig


@memorize()
def sparse_id_pow(x):
    return sparse.identity(2 ** x)


def to_np_if_needed(array):
    return array if isinstance(array, np.ndarray) else array.toarray()

def build_qft_mat_helper(N):
    mat = []
    for i in range(0, N):
        line = [(j * i % N) for j in range(0, N)]
        mat.append(line)
    mat = np.array(mat)
    return mat


def build_qft_mat(n, inv):
    N = 2 ** n
    mat = build_qft_mat_helper(N)
    mul = -1 if inv else 1
    w = np.exp(mul * 2 * np.pi * 1j / N)
    res = (1 / np.sqrt(N)) * (w ** mat)
    return res


def apply_qft_inv(t_a, t_b, state):
    qft_mat_a = build_qft_mat(t_a, inv=True)
    qft_mat = krons(qft_mat_a, np.identity(2 ** t_b))
    res = dots(qft_mat, state)
    return res


def calc_matrix_shift(mat):
    def calc_inner_shift(shift, qubits):
        a = sparse_id_pow(shift)
        b = mat
        if mat.shape == (2, 2):
            mm = 1
        elif mat.shape == (4, 4):
            mm = 2
        else:
            raise Exception(f'wat? expected matrices (2,2) or (4,4), found: {mat.shape}')

        c = sparse_id_pow((qubits - shift - mm))
        return krons(a, b, c)

    return calc_inner_shift


calc_hadamard_shift = calc_matrix_shift(HM_SPARSE)
calc_swap_shift = calc_matrix_shift(SWAP_SPARSE)

hadamard_shift_cache = {}


def hadamard_matrix(shift, qubits):
    key = (shift, qubits)
    mat = hadamard_shift_cache.get(key)
    if mat is None:
        mat = hadamard_shift_cache[key] = calc_hadamard_shift(shift, qubits)
    return mat


def apply_ux(N, x, t_a, t_b, state):
    state = to_np_if_needed(state)
    new_state = np.zeros_like(state)
    for j in range(2 ** t_a):
        old_state_idx = j << t_b
        new_state_idx = old_state_idx + ((x ** j) % N)
        new_state[new_state_idx] = state[old_state_idx]
    return new_state


def build_init_vector(t_):
    state = np.zeros(2 ** t_, dtype=complex)
    state[0] = 1
    return state.reshape(-1, 1)


def check_divider(n, x, r):
    if r % 2 == 1:
        return 0
    xr = x ** (r // 2) % n
    g = gcd(xr + 1, n)
    if 1 < g < n:
        return g
    g = gcd(xr - 1, n)
    if 1 < g < n:
        return g
    return None

def swap_bits_matrix(t_, b1, b2):
    if b1 + 1 != b2:
        raise KeyError("only qubits like (i, i+1) is currently supported!")
    if not (0 <= b1) or not (b2 < t_):
        raise Exception('index out of bound..')
    mat = calc_swap_shift(b1, t_)
    return mat


def phase_gate(m, inv):
    mul = -1 if inv else 1
    phase = np.exp(mul * 2 * np.pi * 1j / (2 ** m))
    return np.array([[1, 0], [0, phase]], dtype=complex)


def block(mat1, mat2):
    assert mat1.shape == (2, 2) and mat2.shape == (2, 2)
    zr = np.zeros_like(mat1)
    crm = np.block([[mat1, zr], [zr, mat2]])
    return crm


# gates_debug = True
gates_debug = False
gates_debug_check = False


class Gate:
    def __init__(self, name, qubits):
        self.name = name
        if qubits is None or isinstance(qubits, list):
            raise Exception()
        self.qubits = qubits

    def apply(self, state):
        raise AttributeError('Not implemented!')

    def collect_lines(self, level):
        return

    def collect_draw_lines(self, level):
        return [''.join(replace_empty(self.draw_idx(i)).rjust(4) for i in range(self.qubits))]

    def draw_idx(self, idx):
        raise AttributeError('Not implemented!')

    def visualize(self, level=1):
        lines = self.collect_draw_lines(level)
        print('\n'.join(lines))


class MatrixGate(Gate):
    def __init__(self, name, qubits, matrix):
        super().__init__(name, qubits)
        self.matrix = matrix

    def apply(self, state):
        res = dots(self.matrix, state)
        return res


class InitZero(Gate):
    def __init__(self, t_):
        super().__init__('Init Zero', t_)
        self.t_ = t_
        self.state = build_init_vector(t_)

    def apply(self, state):
        if state is not None:
            raise KeyError('expected null argument!')
        return self.state

    def draw_idx(self, idx):
        return '?'


def measure_suffix(suffix_size, qubits, state):
    state = to_np_if_needed(state).reshape(-1)
    sfx_mask = (1 << suffix_size) - 1
    probs = np.zeros(2 ** suffix_size)
    for i in range(2 ** qubits):
        # todo optimize: vectorize
        try:
            probs[i & sfx_mask] += abs((state[i]) ** 2)
        except:
            pdb.set_trace()
    print('\tprobs:', len([(i, round(j, 2)) for i, j in enumerate(probs) if not almost_zero(j, 2)]))
    mask = choose_from_probs(probs)
    new_state = np.copy(state)
    for i in range(2 ** qubits):
        if (i & sfx_mask) != mask:
            new_state[i] = 0
    new_state = normalize(new_state)
    new_state = new_state.reshape(-1, 1)
    return mask, new_state


class Measure(Gate):
    def __init__(self, name, qubits, suffix_size):
        super().__init__(name, qubits)
        self.qubits = qubits
        self.suffix_size = suffix_size
        self.name = name

    def apply(self, state):
        selected, res = measure_suffix(self.suffix_size, self.qubits, state)
        self.selected = selected
        return res

    def draw_idx(self, idx):
        th = self.qubits - self.suffix_size
        return 'M' if idx >= th else ''


class Hadamard(MatrixGate):
    def __init__(self, qubits, idx):
        super().__init__('H', qubits, hadamard_matrix(idx, qubits))
        self.qubits = qubits
        self.idx = idx

    def draw_idx(self, idx):
        return '' if idx != self.idx else 'H'


class UxGate_simple(Gate):
    def __init__(self, t_a, t_b, N, x):
        super().__init__('Ux', t_a + t_b)
        self.N = N
        self.x = x
        self.t_a = t_a
        self.t_b = t_b

    def apply(self, state):
        res = apply_ux(self.N, self.x, self.t_a, self.t_b, state)
        return res

    def draw_idx(self, idx):
        return 'x' if idx < self.t_a else 'Ux'


def inf(st):
    return '[' + ', '.join(map(lambda x: str(x)[1:-1], st.reshape(-1).round(4))) + ']'

class Bulk(Gate):
    def __init__(self, name, qubits, gates, post_check=None):
        super().__init__(name, qubits)
        self.gates = gates
        self.post_check = post_check
        self.can_cache_matrix = False
        self.cache_matrix = None

    def matrix(self):
        if self.cache_matrix is None:
            self.cache_matrix = dots(*([g.matrix for g in self.gates][::-1]))
        return self.cache_matrix


    def apply(self, state):
        if self.can_cache_matrix and all(isinstance(g, MatrixGate) for g in self.gates):
            mat = self.matrix()
            # print("using cached!", self.name)
            return dots(mat, state)

        if gates_debug:
            pass
            # print(state.round(2))
        state0 = state
        for g in self.gates:
            if gates_debug:
                print('gate:', g.name, state.shape if state is not None else '???')
                if isinstance(g, MatrixGate):
                    print('matrix:', g.matrix.toarray().round(2))
            # if state is not None and isinstance(state, np.ndarray):
            #     state = csr_matrix(state)
            state = g.apply(state)
            if gates_debug:  # and g.name.startswith('QFT'):
                print('state:')
                print(to_np_if_needed(state).round(4))
                # print(state)
                print()
                # print('\t\t', non_null_amps(state))
            if gates_debug_check:
                if not (isinstance(state, np.ndarray) or isinstance(state, csr_matrix)):
                    raise KeyError(f"Illegal type of state: {type(state)}")
                    # if len(state) != 1:
                    #     raise KeyError('expected vector!')
                    # if not all(isinstance(e, complex) for e in state[0]):
                    #     raise KeyError("all amplitudes should be type complex!")
            if state is None:
                raise KeyError('None is not expected!')
        if gates_debug_check and self.post_check is not None:
            if not self.post_check(state0, state):
                raise KeyError('Post check failed!')
        return state

    def what(self, measure_name):
        measure_gates = [g for g in self.gates if isinstance(g, Measure) and g.name == measure_name]
        if len(measure_gates) != 1:
            raise KeyError(f'1 measure gate with name {measure_name} expected, found: {measure_gates}')
        return measure_gates[0].selected

    def collect_draw_lines(self, level):
        if level > 1:
            return lflatten(g.collect_draw_lines(level - 1) for g in self.gates)
        else:
            return super().collect_draw_lines(level)

    def draw_idx(self, idx):
        return self.name[0]

class MulModGate(Gate):
    def __init__(self, t_a, t_b, y, N, control_qb):
        super().__init__('M_yn', t_a + t_b)
        self.y = y
        self.N = N
        self.t_b = t_b
        self.control_qb = control_qb

    def apply(self, state):
        state = state.toarray().reshape(-1)
        new_state = np.zeros_like(state)
        control_shift = 1 << self.control_qb
        for j in range(2 ** self.t_b):
            pass
        for ii, amp in enumerate(state):
            if (ii & control_shift) != 0:
                new_state_idx = TODO_HERE

    def draw_idx(self, idx):
        if idx == self.control_qb:
            return 'X'
        else:
            return '' if idx < self.t_a else 'M'

class UxGate(Bulk):
    def __init__(self, t_a, t_b, N, x):
        simple_ux = UxGate_simple(t_a, t_b, N, x)
        def post_check(s0, actual):
            expected = simple_ux.apply(s0)
            return expected == actual
        gates = [TODO]
        super().__init__('Ux', t_a + t_b, gates, post_check)

    def apply(self, state):
        return state

class SwapGate(MatrixGate):
    def __init__(self, qubits, qa):
        super().__init__('Swap', qubits, swap_bits_matrix(qubits, qa, qa + 1))
        if not (qa + 1 < qubits):
            raise KeyError(f'illegal swap, expected qa + 1 = {qa + 1} < {qubits} = qubits')
        self.qa = qa

    def draw_idx(self, idx):
        qa = self.qa
        return 'x-' if idx == qa else ('-x' if idx == qa + 1 else '')


def build_inv_phase_matrix(phase_idx):
    phase = np.exp(-2 * np.pi * 1j / (2 ** phase_idx))
    m = np.identity(2, dtype=complex)
    m[1][1] = phase
    return m


@memorize()
def elementary_control_phase(idx):
    phase = build_inv_phase_matrix(idx)
    res = block(E2, phase)
    return csr_matrix(res)


@memorize()
def build_control_phase_near_matrix(qubits, qb, idx):
    m0 = sparse_id_pow(qb)
    m1 = elementary_control_phase(idx)
    m2 = sparse_id_pow((qubits - 2 - qb))
    res = krons(m0, m1, m2)
    return res


class ControlPhaseGateNear(MatrixGate):
    def __init__(self, qubits, qb, idx):
        # TODO rename Inverse,
        super().__init__(f'ControlPhaseNear {idx}', qubits, build_control_phase_near_matrix(qubits, qb, idx))
        self.qb = qb
        self.idx = idx

    def draw_idx(self, idx):
        if idx == self.qb:
            return 'R_' + str(self.idx)
        elif idx == self.qb + 1:
            return 'X'
        else:
            return ''

def build_control_phase(qubits, qb, qb_control, idx):
    if not (qb < qb_control):
        raise KeyError(f'expected qb = {qb} < {qb_control} = qb_control')
    gates = []
    swap_ids = list(range(qb_control - 1, qb - 1, -1))
    #     print('\t\t\t for swap:', qb_control - 1, qb + 1, -1)
    for i in swap_ids:
        gates.append(SwapGate(qubits, i))
    gates.append(ControlPhaseGateNear(qubits, qb, idx))
    for i in swap_ids[::-1]:
        gates.append(SwapGate(qubits, i))
    return gates


class ControlPhaseGate(Bulk):
    def __init__(self, qubits, qb, qb_control, idx):
        super().__init__('ControlPhase' + str(idx), qubits, build_control_phase(qubits, qb, qb_control, idx))
        self.qubits = qubits
        self.qb = qb
        self.qb_control = qb_control
        self.idx = idx
        self.can_cache_matrix = True

    def draw_idx(self, idx):
        if idx == self.qb:
            return 'R_' + str(self.idx)
        elif idx == self.qb_control:
            return 'X'
        else:
            return ''

    def collect_draw_lines(self, level):
        if level == 1:
            return Gate.collect_draw_lines(self, level)
        else:
            return super().collect_draw_lines(level)

class Permutation(MatrixGate):
    def __init__(self, t_a, t_b):
        super().__init__(f'Permutation 2^{t_a}', t_a + t_b, permutation_matrix(t_a, t_b))
        self.t_a = t_a
        self.t_b = t_b

    def draw_idx(self, idx):
        return 'P' if idx < self.t_a else ''


def build_qft_inv_schema(t_a, t_b):
    gates = [Permutation(t_a, t_b)]
    qubits = t_a + t_b
    for cur_qb in range(t_a - 1, -1, -1):
        for cur_other_qb in range(t_a - 1, cur_qb, -1):
            phase_idx = cur_other_qb - cur_qb + 1
            gates.append(ControlPhaseGate(qubits, cur_qb, cur_other_qb, phase_idx))
        gates.append(Hadamard(qubits, cur_qb))
    return gates


def non_null_amps(state: np.ndarray, e=5):
    amps = state.reshape(-1)
    return [i for i in [(i, round(j, e)) for i, j in enumerate(amps)] if i[1] != 0]

class QFT_simple(Gate):
    def __init__(self, t_a, t_b):
        super().__init__('QFT^(-1)', t_a + t_b)
        self.t_a = t_a
        self.t_b = t_b

    def apply(self, state):
        res = apply_qft_inv(self.t_a, self.t_b, state)
        return res

    def draw_idx(self, idx):
        return 'Q' if idx < self.t_a else ''


class QFT(Bulk):
    def __init__(self, t_a, t_b):
        super().__init__('QFT', t_a + t_b, build_qft_inv_schema(t_a, t_b))
        # self.can_cache_matrix = True
        self.t_a = t_a

    def apply(self, state):
        st = super().apply(state)
        return st

    def draw_idx(self, idx):
        return ('Q' if idx < self.t_a else '')


def shor_builder(N, basic_qft: bool, basic_ux: bool):
    def callback(x):
        lgN = log2(N)
        t_a, t_b = int(np.ceil(2 * lgN)), int(np.ceil(lgN))
        print(f'>>N={N}, t_a={t_a}, t_b={t_b}, x={x}, basic_qft={basic_qft}, basic_ux={basic_ux}')
        return ShorConfig(t_a, t_b, N, x, basic_qft, basic_ux)

    return callback


def check_and_return_schema(gate: Bulk):
    if not isinstance(gate, Bulk):
        raise KeyError('Bulk gate expected!')
    if not isinstance(gate.gates[0], InitZero):
        raise KeyError(f'First gate InitZero expected! Found: {gate.gates[0]}')
    if not isinstance(gate.gates[-1], Measure):
        raise KeyError(f'Last gate Measure expected! Found: {gate.gates[-1]}')
    return gate


class Schema:
    def __init__(self, shor_config: ShorConfig, schema: Bulk):
        self.shor_config = shor_config
        self.schema = check_and_return_schema(schema)

    def run(self):
        self.schema.apply(None)

    def what(self, measure_key: str):
        return self.schema.what(measure_key)

    def visualize(self, level=2):
        inc = 2 if self.shor_config.t_ < 50 else 1
        self.schema.visualize(level)


def build_schema(sc: ShorConfig):
    def hadamars_post_check(s0, actual):
        actual = to_np_if_needed(actual)
        return len(non_null_amps(actual)) == 2 ** sc.t_a
    qft_constructor = QFT_simple if sc.basic_qft else QFT
    ux_constructor = UxGate_simple if sc.basic_ux else UxGate
    schema = Bulk('Shor', sc.t_, [
        InitZero(sc.t_),
        Bulk('Hadamars', sc.t_, [Hadamard(sc.t_, i) for i in range(0, sc.t_a)], hadamars_post_check),
        ux_constructor(sc.t_a, sc.t_b, sc.N, sc.x),
        Measure('t_b', sc.t_, sc.t_b),
        qft_constructor(sc.t_a, sc.t_b),
        Measure('t_', sc.t_, sc.t_)
    ])
    return Schema(sc, schema)


@timeit
def try_find_0(sc: ShorConfig):
    schema = build_schema(sc)
    schema.run()
    selected_b, selected = schema.what('t_b'), schema.what('t_')
    if (selected & ((1 << sc.t_b) - 1)) != selected_b:
        msg = f'wat???, {bin(selected)}, {bin(selected & ((1 << sc.t_b) - 1))}, {bin(selected_b)}'
        raise KeyError(msg)
    ch = selected >> sc.t_b
    r = Fraction(ch, 1 << sc.t_a).denominator
    res = check_divider(sc.N, sc.x, r)
    return None if res == 0 else res

def do_attempts(shor_config, attempts):
    for at in range(attempts):
        res = try_find_0(shor_config)
        msg = 'FAIL' if res is None else 'SUCCESS!!' + '='*20 + f'> {res}'
        print(f'attempt {at+1}... ', msg)
        if res is not None:
            return True
    return False


def try_find(N, attempts, basic_qft, basic_ux):
    sc_without_x = shor_builder(N, basic_qft, basic_ux)
    for x in range(N - 1, 1, -1):
        if gcd(x, N) != 1:
            continue
        sc = sc_without_x(x)
        print('x:', x)
        success = do_attempts(sc, attempts)
        if success:
            return
    print('FAIL...')


def try_find_2(N, x, basic_qft, basic_ux):
    if gcd(x, N) != 1:
        raise KeyError(f'expected that gcd(N={N}, x={x})=1!')
    sc = shor_builder(N, basic_qft, basic_ux)(x)
    success = do_attempts(sc, 15)
    if not success:
        print('FAIL...')

class Helper:
    def __init__(self, N=15, x=None, basic_qft=False, basic_ux=True, debug=False, debug_check=False):
        self.N = N
        self.x = x
        self.basic_qft = basic_qft
        self.basic_ux = basic_ux
        global gates_debug, gates_debug_check
        gates_debug = debug
        gates_debug_check = debug_check
        print('debug:', debug)
        print('debug check:', debug_check)
        

    def run(self, attempts=10):
        if self.x is None:
            try_find(self.N, attempts, self.basic_qft, self.basic_ux)
        else:
            try_find_2(self.N, self.x, self.basic_qft, self.basic_ux)

    def draw(self, lev=1):
        sc = shor_builder(self.N, self.basic_qft, self.basic_ux)(self.x)
        build_schema(sc).visualize(lev)

if __name__ == '__main__':
    fire.Fire(Helper)
    # print(build_qft_mat(6, True))
