# shor_algorithm_simulator
Shor algorithm simulator in Python

## Install
```bash
pip install -r requirements.txt
```

## Run examples
Find the divisors of 15
```bash
python main.py run --N=15
```

Visualize quantum scheme for 35
```bash
python main.py draw --N=35 --lev=3
```

| lev | details                                            |
|-----+----------------------------------------------------|
|   1 | Nothing                                            |
|   2 | High-level gates (Hadamard, Ux, Measure, QFT^(-1)) |
|   3 | Low-level gates (expanded QFT^(-1))                |
|   4 | Lowest-level gates (expanded Control-Phase gates)  |
