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

g = dpalign.align_and_display(seq1, seq2, seq1 + "." + seq2 + ".png", False)

print("Image saved to " + seq1 + "." + seq2 + ".png")


narration = "This animation shows the alignment of sequences " + " ".join(
     [a for a in seq1]) + " and " + " ".join([a for a in seq2]) + ". "
narration += "Dynamic programming is used to calculate the best route from " \
             "the origin to each node. The arrows " \
             "indicate the source of the optimal path for that node. " \
             "Optimal alignments are found by following those traceback " \
             "arrows from the final sink back to the origin. "

narration += ("In this particular case, " +
              str(len(dpalign.map_alignments(g, seq1, seq2))) +
              " alignments are equally good.")

dpalign.createAudioFile("output/narration", narration)

dpalign.saveMovie(seq1 + "." + seq2 + ".png.",
                  "theoutput", 3, "output/narration.aiff")
