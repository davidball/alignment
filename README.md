# alignment

Implementation of DP Recursion graphic generation using networkx.

Requirements:
Python Version 3.5
MacOS for the "say" command is required for the audio narration.


Sample Usage:
	python alignment_runner.py GWWPDT WRRKHY

Any two protein sequences can be supplied in place of the command line arguments.  If none are supplied it defaults to using those two sequences as an example. 

Output
Creates an image file named with the provided sequences in format SEQ1.SEQ2.png
Creates a narrated video demonstrating the key steps in the alignment algorithm with file named with the provided sequences in format SEQ1.SEQ2.mp4
Outputs the best alignments to the console. 

[Sample video output](https://www.youtube.com/watch?v=QkxCH49Ry2E)

See also the jupyter notebook alignment_graph.ipynb for examples of dynamically changing the scoring function.