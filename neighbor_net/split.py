
__author__ ="Daniel H. Huson"

class Split:
    def __init__(self, partA, partB, weight=1):
        self.partA=list(partA)
        self.partB = list(partB)
        self.weight=weight

    def __str__(self):
        return f'{self.partA} {self.partB} {self.weight}'

    def partA(self):
        return self.partA

    def partB(self):
        return self.partB

    def size(self):
        return min(len(self.partA), len(self.partB))

    def trivial(self):
        return self.size()==1

    def weight(self):
        return self.weight

    def set_weight(self,weight):
        self.weight=weight