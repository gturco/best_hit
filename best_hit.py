from bed_utils import BlastLine
import itertools
from collections import defaultdict
from flatfeature import Bed

def group_hits(s):
    """groups values to same key in a list"""
    d = defaultdict(list)
    for k, v in s:
        d[k].append(v)
    return d


def blast_grouped(blast_file,query):
    """groups all hits to key sorts by largest hit and returns best hit
    and score """
    fp = file(blast_file)
    blasts = sorted([BlastLine(line) for line in fp], \
                             key=lambda b: b.score, reverse=True)
    if query: blast_list = [(b.query, (b.score,b.subject)) for b in blasts if b.evalue < 1e-10]
    else: blast_list = [(b.subject, (b.score,b.query)) for b in blasts if b.evalue < 1e-10]
    blast_grouped = group_hits(blast_list)
    return blast_grouped

def write_best_hits(best_hits, out_file):
    write_file = open(out_file,"wb")
    for k in best_hits:
        hits = best_hits[k]
        hits.sort()
        score,subject = hits[-1]
        w = "{0}\t{1}\t{2}\n".format(k,subject,score)
        write_file.write(w)


#best_hits = best_hit("/Users/gt/data/paper4/rice_j_setaria_n/rice_j_setaria_n.blast")
#write_best_hits(best_hits,"/Users/gt/data/paper4/rice_j_setaria_n/rice_j_setaria_n.best_hit",Bed("/Users/gt/setaria_n.bed"))
#best_hits = best_hit("/Users/gt/data/paper4/rice_j_sorghum_n/rice_j_sorghum_n.blast")
#write_best_hits(best_hits,"/Users/gt/data/paper4/rice_j_sorghum_n/rice_j_sorghum_n.best_hit",Bed("/Users/gt/sorghum_n.bed"))

best_hits = blast_grouped("/Users/gt/data/paper4/rice_j_sorghum_n/rice_j_sorghum_n.blast",True)
write_best_hits(best_hits,"/Users/gt/rice_j_sorghum_n.best_hit")



