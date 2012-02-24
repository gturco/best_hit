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

def best_hit(blast_file):
    """groups all hits to key sorts by largest hit and returns best hit and
    score """
    fp = file(blast_file)
    #blast = BlastLine(blast_file)
    blasts = sorted([BlastLine(line) for line in fp], \
                        key=lambda b: b.score, reverse=True)

    blast_list = [(b.subject, (b.score,b.query)) for b in blasts if b.evalue < 1e-10]
    blast_grouped = group_hits(blast_list)
    #blast_grouped = itertools.groupby(blast_list, key=lambda x:x[0])
    best_hits = dict((k,max(blast_grouped[k])) for k in blast_grouped)
    return best_hits

def write_best_hits(best_hits, out_file, sbed):
    write_file = open(out_file,"wb")
    for k in best_hits:
        score,subject = best_hits[k]
        sfeat = sbed.accn(subject)
        w = "{0}\t{1}\t{2}\t{3}\t{4}\n".format(k,subject,sfeat["start"],sfeat["end"],score)
        write_file.write(w)


#best_hits = best_hit("/Users/gt/data/paper4/rice_j_setaria_n/rice_j_setaria_n.blast")
#write_best_hits(best_hits,"/Users/gt/data/paper4/rice_j_setaria_n/rice_j_setaria_n.best_hit",Bed("/Users/gt/setaria_n.bed"))
#best_hits = best_hit("/Users/gt/data/paper4/rice_j_sorghum_n/rice_j_sorghum_n.blast")
#write_best_hits(best_hits,"/Users/gt/data/paper4/rice_j_sorghum_n/rice_j_sorghum_n.best_hit",Bed("/Users/gt/sorghum_n.bed"))

best_hits = best_hit("/Users/gt/data/paper4/sorghum_n_setaria_n.blast")
write_best_hits(best_hits,"/Users/gt/data/paper4/setaria_n.best_hit",Bed("/Users/gt/sorghum_n.bed"))



