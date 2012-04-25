from bed_utils import BlastLine,Bed
from best_hit_ortho import get_pairs,get_pos

"""checks for syn to left and right of gene"""


def main(qbed_file,sbed_file,pairs_file,out,query=True):
    out_fh = open(out,"wb")
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
    for qi in qaccns:
        ### had to change from int search to term because of issues with
        ### merging see Os12g12370
        left_ortho = qbed[qi-1].accn in qorthos
        right_ortho = qbed[qi+1].accn in qorthos
        line = "{0}\t{1}\t{2}\n".format(qbed[qi].accn,left_ortho,right_ortho)
        out_fh.write(line)
    out_fh.close()

#main(qbed_nolocaldups,sbed_nolocaldups,pair_file,true/flase)
main("/Users/gt/new/rice_j_set/rice_j.nolocaldups.with_new.local","/Users/gt/new/rice_j_set/setaria_n.nolocaldups.with_new.local","/Users/gt/new/rice_j_set/rice_j_setaria_n.pairs.txt.local","/Users/gt/new/rice_j_set/rice_j_syn.txt",True)
#main("/Users/gt/new/rice_j.nolocaldups.with_new.local","/Users/gt/new/sorghum_n.nolocaldups.with_new.local","/Users/gt/new/rice_j_sorghum_n.pairs.txt.local","/Users/gt/new/rice_j_syn.txt",True)

