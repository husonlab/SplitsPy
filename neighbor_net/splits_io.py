import sys
import csplit

__author__ = "Daniel H. Huson"


def print_splits_nexus(labels: [str], splits: [csplit.CSplit], cycle: [int], filename="-") -> None:
    if filename == "-":
        outs = sys.stdout
    else:
        outs = open(filename, mode="w")

    print("#nexus",file=outs)

    print("BEGIN taxa;",file=outs)
    print("DIMENSIONS nTax=",len(labels),";",sep="",file=outs)
    print("TAXLABELS",file=outs)
    for label in labels:
        print("'",label,"'",sep="", file=outs)
    print(";",file=outs)
    print("END;",file=outs)

    print("BEGIN SPLITS;",file=outs)
    print("DIMENSIONS nTax=",len(labels)," nSplits=",len(splits),";",sep="",file=outs)
    print("FORMAT labels=no weights=yes confidences=no;",file=outs)
    # print("PROPERTIES...)
    print("CYCLE",end="",file=outs)
    for i in range(1,len(cycle)):
        print(" ",cycle[i], sep="",end="", file=outs)
    print(";", file=outs)
    print("MATRIX",file=outs)
    for split in splits:
        print(f'{split.get_weight():.8f}', end="\t", file=outs)
        first = True
        for t in split.part1():
            if first:
                first = False
            else:
                print(" ",end="",file=outs)
            print(t,end="", file=outs)
        print(",", file=outs)
    print(";",file=outs)
    print("END;",file=outs)

    if outs != sys.stdout:
        outs.close()










