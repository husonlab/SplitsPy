__author__ = "Daniel H. Huson"


class Split:
    def __init__(self, n_tax: int, part_a: [int], weight=1):
        self.part_a = set(part_a)
        self.part_b = set()
        for t in range(1, n_tax + 1):
            if t not in self.part_a:
                self.part_b.add(t)
        self.weight = weight

    def __str__(self):
        return f'{self.part_a} {self.part_b} {self.weight}'

    def get_part_a(self) -> [int]:
        return self.part_a

    def get_part_b(self) -> [int]:
        return self.part_b

    def get_1(self) -> [int]:
        return self.part_a if 1 in self.part_a else self.part_b

    def get_size(self) -> int:
        return min(len(self.part_a), len(self.part_b))

    def is_trivial(self) -> bool:
        return self.get_size() == 1

    def get_weight(self) -> float:
        return self.weight

    def set_weight(self, weight: float) -> None:
        self.weight = weight
