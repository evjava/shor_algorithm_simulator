# shor_algorithm_simulator
Shor algorithm ( https://en.wikipedia.org/wiki/Shor%27s_algorithm ) simulator in Python

Circuit is based on the materials of "6.3 Shor’s Factoring Algorithm for Breaking RSA Quantumly" from book "Explorations in Quantum Computing" ( https://link.springer.com/book/10.1007/978-1-84628-887-6 )

Note that in this implementation UX gate is not decomposed to the basic quantum gates.

## Install
```bash
pip install -r requirements.txt
```

## Run examples

### Run shor algorithm for N=15
```bash
python main.py run --N=15
```

### Visualize quantum circuit for N=35
```bash
python main.py draw --N=35 --lev=3
```

Other lev options:
| lev | details                                            |
| --- | -------------------------------------------------- |
|   1 | Nothing                                            |
|   2 | High-level gates (Hadamard, Ux, Measure, QFT^(-1)) |
|   3 | Low-level gates (expanded QFT^(-1))                |
|   4 | Lowest-level gates (expanded Control-Phase gates)  |

### Run Shor algorithm for N=15, use ready matrix for QFT^(-1)
```bash
python main.py run --N=35 --lev=3 --basic-qft=True
```

### Run Shor algorithm for N=21, use x=5
```bash
python main.py run --N=21 --x=5
```
