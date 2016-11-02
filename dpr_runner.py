import sys
import dpalign

if len(sys.argv) >= 3:
    seq1 = sys.argv[1]
    seq2 = sys.argv[2]
else:
    seq1 = 'GWWPDT'
    seq2 = 'WRRKHY'
    print("No sequence data provided at the command line, using the defaults")


print("Aligning %s and %s" % (seq1, seq2))

show_plot = False
generate_movie = True
g = dpalign.align_and_display(seq1, seq2, seq1 + "." + seq2 + ".png",
                              show_plot, generate_movie)

print("Image saved to " + seq1 + "." + seq2 + ".png")

