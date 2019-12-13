import os
import sys
from RNASeqReader import RNAseqReader
from shared import Interval
RNAseqReader = RNAseqReader(fname="Z:/scratch/202001051010polina_data/3DPredictor/input/K562/RNA-seq/rna-seqPolyA.tsvpre.txt",
                                           name="RNA")
RNAseqReader.read_file(rename={ "Gene name": "gene",
                              "Gene start (bp)": "start",
                              "Gene end (bp)": "end",
                              "Chromosome/scaffold name": "chr",
                              "FPKM": "sigVal"},
                      sep="\t")
interval = Interval("chr14", 19110200, 19253800)
print(RNAseqReader.chr_data["chr14"])
interval_data = RNAseqReader.get_interval(interval)
print(interval_data)

# RNAseqPG = SmallChipSeqPredictorGenerator(RNAseqReader,
#                                                   window_size=params.window_size,
#                                                   N_closest=3)
