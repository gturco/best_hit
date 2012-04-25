from bed_utils import BlastLine,Bed
from  itertools import product
from collections import defaultdict
from best_hit import blast_grouped


def get_same_seqid(left_ortho,right_ortho,hits,sbed):
    """ if hit same seqid and falls inbetween orthos then get best hit score"""
    ortho_hits = []
    seqid = left_ortho.seqid
    for score,subject in hits:
        hit = sbed[subject][1]
        if hit.seqid == seqid:
            if hit.start > left_ortho.start and hit.end < right_ortho.end:
                ortho_hits.append((score,subject))
    ortho_hits.sort()
    try:
        return ortho_hits[-1][1]
    except IndexError:
        return None

def get_diff_seqid(left_ortho,right_ortho,hits,sbed):
    """ if hits are not the same seqid then it finds best hit on either
    seqid near region"""
    ortho_hits = []
    for score,subject in hits:
        hit = sbed[subject][1]
        if hit.seqid == left_ortho.seqid and hit.start > left_ortho.start:
                ortho_hits.append((score,subject))
        elif hit.seqid == right_ortho.seqid and hit.end < right_ortho.end:
            ortho_hits.append((score,subject))
    ortho_hits.sort()
    try:
        return ortho_hits[-1][1]
    except IndexError:
        return None

def get_best_hits(left_qaccn,right_qaccn,blast,qbed,sbed,pairs,padding):
    """finds all genes in between ortho region and finds their best ortho hit"""
    genes = []
    for iqaccn in range(left_qaccn+1,right_qaccn):
        qaccn = qbed[iqaccn]
        if qaccn.accn not in blast.keys(): continue
        hits = blast[qaccn.accn]
        left_ortho_list = pairs[qbed[left_qaccn].accn]
        right_ortho_list = pairs[qbed[right_qaccn].accn]
        for left_ortho,right_ortho in product(left_ortho_list,right_ortho_list):
            left_ortho = sbed[left_ortho][1]
            right_ortho = sbed[right_ortho][1]
            left_ortho.start = max(0,left_ortho.start - padding)
            right_ortho.end = right_ortho.end + padding

            if left_ortho.seqid == right_ortho.seqid:
                best_hit = get_same_seqid(left_ortho,right_ortho,hits,sbed)
                genes.append((qaccn.accn,best_hit))
            else:
                best_hit = get_diff_seqid(left_ortho,right_ortho,hits,sbed)
                genes.append((qaccn.accn,best_hit))
    return genes

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
    for gene,hit in best_hits:
        w = "{0}\t{1}\n".format(gene,hit)
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
    blasts = blast_grouped(blast_file,query)
    qaccns.sort()
    best_hits = []
    for qi,q in enumerate(qaccns):
        if qi == 0: continue
        if qi == len(qaccns): continue
        left_ortho = qaccns[qi-1]
        right_ortho = qaccns[qi]
        if (right_ortho - left_ortho) == 1: continue
        #### if new ortho - old is == 1 no orthos inbetween
        new_hits = get_best_hits(left_ortho,right_ortho,blasts,qbed,sbed,pairs_dict,padding)
        best_hits += new_hits


    write_best_hits(out_fh,best_hits)
    return best_hits


main("/Users/gt/data/paper4/rice_j.bed","/Users/gt/data/paper4/setaria_n.bed","/Users/gt/data/paper4/rice_j_setaria_n/rice_j_setaria_n.blast","/Users/gt/data/paper4/rice_j_setaria_n/rice_j_setaria_n.pairs.txt","/Users/gt/data/paper4/rice_j_setaria_n/rice_j_setaria_n.pairs.30000",30000)
#best_hits = best_hit("/Users/gt/dick_m_tair_10.blast")
#write_best_hits(best_hits,"/Users/gt/dick_m_tair_10.best_hits",Bed("/Users/gt/tair_10.bed"))
