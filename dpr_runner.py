import sys
import dpalign
import networkx as nx
import matplotlib.pyplot as plt

if len(sys.argv)>=3:
    seq1 = sys.argv[1]
    seq2 = sys.argv[2]
else:
    seq1='GWWPDT'
    seq2='WRRKHY'
    print("No sequence data provided at the command line, so using the defaults")


print("Aligning %s and %s" % (seq1,seq2))

dpalign.align_and_display(seq1,seq2,seq1 + "." + seq2 + ".png",False)

print("Image saved to " + seq1 +"."+seq2+".png")