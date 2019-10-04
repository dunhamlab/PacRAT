# this file aligns barcodes

from Bio.Align.Applications import MuscleCommandline
from StringIO import StringIO
from Bio import AlignIO

muscle_exe = r"/net/gs/vol3/software/modules-repo/RHEL6/muscle"

inputfile = 
outputfile = "test_alignment.fasta"