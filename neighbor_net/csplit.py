__author__ = "Daniel H. Huson"


class CSplit:
    def __init__(self, n_tax: int, cycle: [int], start_pos: [int], end_pos: int,  weight=1):
        self.n_tax = n_tax
        self.cycle = cycle
        self.start_pos = start_pos
        self.end_pos = end_pos
        self.weight = weight

    def __str__(self):
        return f'{self.cycle[self.start_pos]} - {self.cycle[self.end_pos]} {self.weight}'

    def get_start_pos(self) -> int:
        return self.start_pos

    def get_end_pos(self) -> int:
        return self.end_pos

    def get_1(self) -> [int]:
        a = [self.cycle[i] for i in range(1,self.start_pos)]
        a.extend(self.cycle[i] for i in range(self.end_pos+1, self.n_tax+1))
        a.sort()
        return a

    def get_2(self) -> [int]:
        a = [self.cycle[i] for i in range(self.start_pos,self.end_pos+1)]
        a.sort()
        return a

    def get_size(self) -> int:
        diff = self.end_pos-self.start_pos+1
        return min(diff, self.n_tax-diff)

    def is_trivial(self) -> bool:
        return self.get_size() == 1

    def get_weight(self) -> float:
        return self.weight

    def set_weight(self, weight: float) -> None:
        self.weight = weight
