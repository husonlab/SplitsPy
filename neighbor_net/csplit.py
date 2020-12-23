__author__ = "Daniel H. Huson"


class CSplit:
    def __init__(self, n_tax: int, cycle: [int], start_pos: [int], end_pos: int,  weight=1):
        self.__n_tax = n_tax
        self.__cycle = cycle

        if start_pos == 1:
            start_pos = end_pos+1
            end_pos = n_tax

        self.__start_pos = start_pos
        self.__end_pos = end_pos
        self.weight = weight

    def __str__(self):
        return f'{self.__start_pos} - {self.__end_pos} ({self.__cycle[self.__start_pos]} - {self.__cycle[self.__end_pos]}) {self.weight: .8f}'

    def start_pos(self) -> int:
        return self.__start_pos

    def end_pos(self) -> int:
        return self.__end_pos

    def part1(self) -> [int]:
        a = [self.__cycle[i] for i in range(1, self.__start_pos)]
        a.extend(self.__cycle[i] for i in range(self.__end_pos + 1, self.__n_tax + 1))
        a.sort()
        return a

    def part2(self) -> [int]:
        a = [self.__cycle[i] for i in range(self.__start_pos, self.__end_pos + 1)]
        a.sort()
        return a

    def size(self) -> int:
        diff = self.__end_pos - self.__start_pos + 1
        return min(diff, self.__n_tax - diff)

    def is_trivial(self) -> bool:
        return self.size() == 1

    def get_weight(self) -> float:
        return self.weight

    def set_weight(self, weight: float) -> None:
        self.weight = weight
