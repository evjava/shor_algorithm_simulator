import fire

from main import *


class QFTRunner:
    def __init__(self):
        self.is_small = True

    def draw(self, N=15, lev=2):
        if self.is_small:
            schema = QFT(N, 0)
        else:
            schema = QFT_simple(N, 0)
        schema.visualize(lev)


if __name__ == '__main__':
    fire.Fire(QFTRunner)
