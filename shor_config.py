from dataclasses import dataclass


@dataclass
class ShorConfig:
    t_a: int
    t_b: int
    N: int
    x: int
    basic_qft: bool
    basic_ux: bool

    @property
    def t_(self):
        return self.t_a + self.t_b
