#!/usr/bin/python 
# Universitat Potsdam
# Author Gaurav Sablok
# date: 2024-2-29
def gfafastawrite(gfafile, filewrite):
    """
    a gfa to fasta write for writing the graph connections for the 
    GFA files
    """
    sort = list(pd.DataFrame([line.split("\t") for line in open("/home/gaurav/Downloads/MT.gfa")], \
                                columns = ["a","b","c","d","e","f","g"])["c"])
    indices = [i for i in range(len(sort)) if "+" not in sort[i] and "-" not in sort[i]]
    sequences = [i for i in sort if "+" not in i and "-" not in i]
    ids = list(pd.DataFrame([line.split("\t") for line in open("/home/gaurav/Downloads/MT.gfa")], \
                            columns = ["a","b","c","d","e","f","g"])["b"])[0:len(indices)]
    with open(filewrite, "w") as fastawrite:
        for i in range(len(ids)):
            fastawrite.write(">"+ids[i])
            fastawrite.write(sequences[i])
    fastawrite.close()
