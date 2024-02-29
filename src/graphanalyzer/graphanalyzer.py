import pandas as pd
import searborn as sns
sns.set_theme()


def pafplotter(paffile = None, \
                    querychr = None, \
                        mapping = None, \
                                lengthplot = "FALSE", \
                                        catplot = "FALSE", \
                                                    mapping_quality_plot = "FALSE", \
                                                                scatterplot = "FALSE"):
    # Universitat Potsdam
    # Author Gaurav Sablok
    # date: 2024-2-28
    """"
    a pafplotter function which will allow you to plot many analysis
    and summary stats from the paf alignment files for the graphs. It 
    analyses the paf alignment files faster and you can chunck the indiviual
    block of the code as the separate functions. 
    """
    if paffile is not None and querychr is None:
        pafdataframe = pd.DataFrame([line.strip().split("\t")[0:12] for line in open(paffile)], 
        columns=["query", "query_length", "query_start", "query_end", "strand", "target", \
                                 "target_length", "target_start", "target_end", "residue_matches", 
                                                                "alignment_length", "mapping_quality"])
        return pafdataframe
    if paffile and querychr is not None:
        pafdataframe = pd.DataFrame([line.strip().split("\t")[0:12] for line in open(paffile)], 
        columns=["query", "query_length", "query_start", "query_end", "strand", "target", \
                                 "target_length", "target_start", "target_end", "residue_matches", 
                                                                "alignment_length", "mapping_quality"])
        selectedpafdataframe = pafdataframe.where(pafdataframe["query"] == querychr).dropna()
        return selectedpafdataframe
    if paffile and mapping is not None:
        pafdataframe = pd.DataFrame([line.strip().split("\t")[0:12] for line in open(paffile)], 
        columns=["query", "query_length", "query_start", "query_end", "strand", "target", \
                                 "target_length", "target_start", "target_end", "residue_matches", 
                                                                "alignment_length", "mapping_quality"])
        selectedpafdataframe = pafdataframe.where(pafdataframe["mapping_quality"] == mapping).dropna()
        return selectedpafdataframe
    if paffile is not None and querychr is not None and lengthplot == "TRUE":
        pafdataframe = pd.DataFrame([line.strip().split("\t")[0:12] for line in open(paffile)], 
        columns=["query", "query_length", "query_start", "query_end", "strand", "target", \
                                 "target_length", "target_start", "target_end", "residue_matches", 
                                                                "alignment_length", "mapping_quality"])
        selectedpafdataframe = pafdataframe.where(pafdataframe["query"] == querychr).dropna()
        pd.DataFrame(pafdataframe.where(pafdataframe["query"] == querychr).dropna()["query_end"].apply(lambda num: int(num)) - \
               pafdataframe.where(pafdataframe["query"] == querychr).dropna()["query_start"].apply(lambda num: int(num)), \
                                                                                    columns = ["plot"]).plot(kind = "bar")
        return selectedpafdataframe
    if paffile is not None and querychr is not None and catplot == "TRUE":
        pafdataframe = pd.DataFrame([line.strip().split("\t")[0:12] for line in open(paffile)], 
        columns=["query", "query_length", "query_start", "query_end", "strand", "target", \
                                 "target_length", "target_start", "target_end", "residue_matches", 
                                                                "alignment_length", "mapping_quality"])
        selectedpafdataframe = pafdataframe.where(pafdataframe["query"] == querychr).dropna()
        catcorrdinates = pd.DataFrame(pafdataframe.where(pafdataframe["query"] == querychr).dropna()["query_end"].apply(lambda num: int(num)) - \
               pafdataframe.where(pafdataframe["query"] == querychr).dropna()["query_start"].apply(lambda num: int(num)), \
                                                                                    columns = ["plot"])
        return sns.catplot(catcorrdinates)
    if paffile is not None and mapping_quality_plot == "TRUE":
        pafdataframe = pd.DataFrame([line.strip().split("\t")[0:12] for line in open(paffile)], 
        columns=["query", "query_length", "query_start", "query_end", "strand", "target", \
                                 "target_length", "target_start", "target_end", "residue_matches", 
                                                                "alignment_length", "mapping_quality"])
        return sns.catplot(pafdataframe, kind = "box", x="residue_matches", y="alignment_length", height = 10)
    if paffile is not None and scatterplot == "TRUE":
        pafdataframe = pd.DataFrame([line.strip().split("\t")[0:12] for line in open(paffile)], 
        columns=["query", "query_length", "query_start", "query_end", "strand", "target", \
                                 "target_length", "target_start", "target_end", "residue_matches", 
                                                                "alignment_length", "mapping_quality"])
        return sns.scatterplot(pafdataframe, x="query", y="mapping_quality")
    

def locationfetch(paffile,querychr):
    """"
    a location fetcher for the paf alignment as a 
    list with the coordinates. Here you can provide
    the query of the alignment file.
    """

    if paffile is not None and querychr is not None:
        paf = paffile
        query = querychr
        pafdataframe = pd.DataFrame([line.strip().split("\t")[0:12] for line in open(paf)], 
        columns=["query", "query_length", "query_start", "query_end", "strand", "target", \
                                 "target_length", "target_start", "target_end", "residue_matches", 
                                                                "alignment_length", "mapping_quality"])
        querycorrdinates = list(pafdataframe.where(pafdataframe["query"] == querychr).dropna().iloc())
        return querycorrdinates
    
def filterstrand(paffile, strand):
    """"
    fetching the alignment length on the query strand 
    for the specific strand. It takes the stand as the 
    str and then locates across the dataframe and then
    map the int to the str conversion. 
    """
    if paffile is not None and strand is not None:
        paf = paffile
        querystrand = str(strand)
        pafdataframe = pd.DataFrame([line.strip().split("\t")[0:12] for line in open(paf)], \
                    columns=["query", "query_length", "query_start", "query_end", "strand", "target", \
                                        "target_length", "target_start", "target_end", "residue_matches", 
                                                                         "alignment_length", "mapping_quality"])
        start = list(map(lambda num: int(num),list(pafdataframe. \
                                                         where(pafdataframe["strand"] == "+").dropna()["query_start"])))
        end = list(map(lambda num: int(num),list(pafdataframe. \
                                                      where(pafdataframe["strand"] == "+").dropna()["query_end"])))
        querystrandalign = []
        for i in range(len(start)):
            querystrandalign.append(end[i] - start[i])
        querylength = pd.DataFrame(querystrandalign, columns = ["alignedlength"]).plot(kind = "bar")
    return querylength
