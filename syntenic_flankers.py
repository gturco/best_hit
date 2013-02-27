from bed_utils import BlastLine,Bed
from  itertools import product
from collections import defaultdict
from best_hit import blast_grouped


def get_pos(qorder,qorthos):
    """returns a list of numbers of genes with orthologs"""
    return [qorder[accn][0] for accn in qorthos]

def get_pairs(pair_file,query):
    pairs = []
    pairs_dict = defaultdict(list)
    for line in open(pair_file):
        if line[0] == "#": continue
        line = line.strip().split("\t")
        qaccn,saccn = line
        pairs.append((qaccn,saccn))
        if query: pairs_dict[qaccn].append(saccn)
        else: pairs_dict[saccn].append(qaccn)
    return pairs,pairs_dict

def write_best_hits(out_fh,best_hits):
    write_file = open(out_fh,"wb")
    for gene in best_hits:
        w = "{0}\n".format(gene)
        write_file.write(w)
    write_file.close()

def main(qbed_file,sbed_file,blast_file,pairs_file,out_fh,padding,query=True):
    if query:
        qbed = Bed(qbed_file)
        sbed = Bed(sbed_file).get_order()
    else:
        qbed = Bed(sbed_file)
        sbed = Bed(qbed_file).get_order()

    qorder = qbed.get_order()
    pairs,pairs_dict = get_pairs(pairs_file,query)
    if query: qorthos = list(set([qaccn for qaccn,saccn in pairs]))
    else: qorthos = list(set([saccn for qaccn,saccn in pairs]))
    qaccns = get_pos(qorder,qorthos)
    qaccns.sort()
    flankers = []
    for qi,q in enumerate(qaccns):
        if qi == 0: continue
        if qi == len(qaccns) - 1: continue
        left_pos = q -1
        right_pos = q + 1
        if left_pos in qaccns and right_pos in qaccns:
            #print qbed[q].accn
            flankers.append(qbed[q].accn)
    write_best_hits(out_fh,flankers)
    return flankers

#### run with sorghum .... query = flase
######### in correct FIX
main("/Users/gt/new/rice_j.nolocaldups.bed.with_new.all.local","/Users/gt/new/sorghum_nn.nolocaldups.bed.with_new.all.local","/Users/gt/new/rice_j_sorghum_nn.blast","/Users/gt/new/rice_j_sorghum_nn.pairs.txt.local","/Users/gt/new/rice_j_sorghum.flankers2",30000)
#best_hits = best_hit("/Users/gt/dick_m_tair_10.blast")
#write_best_hits(best_hits,"/Users/gt/dick_m_tair_10.best_hits",Bed("/Users/gt/tair_10.bed"))
