#!/usr/bin/python 
# Universitat Potsdam
# Author Gaurav Sablok
# date: 2024-2-29
def gfaformatter(gfafile):
    import pandas as pd
    import seaborn as sns
    import matplotlib as plt
    import re
    if gfafile is not None:
        gfa = gfafile
        # add the array of the columnames
        columnnames = []
        readgfa = pd.read_csv([line.split("\t") for line in open(gfa)], columns = columnsnames)
        sequencematch = []
        for i in readgfa.iterrows():
            for j in range(len(["A", "T", "G", "C"])):
                if j in i:
                    sequencematch.append(i)
        print(f"the number of the sequence elements in the gfa are: len(sequencematch)")
