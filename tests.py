# import os
# import sys
# from RNASeqReader import RNAseqReader
# from shared import Interval
# import pickle
# RNAseqReader = RNAseqReader(fname="Z:/scratch/202001051010polina_data/3DPredictor/input/K562/RNA-seq/rna-seqPolyA.tsvpre.txt",
#                                            name="RNA")
# RNAseqReader.read_file(rename={ "Gene name": "gene",
#                               "Gene start (bp)": "start",
#                               "Gene end (bp)": "end",
#                               "Chromosome/scaffold name": "chr",
#                               "FPKM": "sigVal"},
#                       sep="\t")
# interval = Interval("chr14", 19110200, 19253800)
# print(RNAseqReader.chr_data["chr14"])
# interval_data = RNAseqReader.get_interval(interval)
# print(interval_data)

# RNAseqPG = SmallChipSeqPredictorGenerator(RNAseqReader,
#                                                   window_size=params.window_size,
#                                                   N_closest=3)
# data = pickle.load(open("Z:/scratch/202001051010polina_data/3DPredictor/input/sequence/hg38/chr22.fa.chr22.fa['A', 'C', 'G', 'N', 'T', 'a', 'c', 'g', 'n', 't']default_chr_name_ranamer","rb"))
# print(data.chrmSizes)
import straw
# result = straw.straw("KR", "Z:/scratch/202001051010polina_data/3DPredictor/input/H1/4DNFI2TK7L2F.hic", "21", "21", "BP", 2500000)
# result = straw.straw("KR", "https://hicfiles.s3.amazonaws.com/hiseq/gm12878/in-situ/combined.hic", "14", "14", "BP", 2500000)
# print(result)
result = straw.straw("KR", "/mnt/scratch/ws/psbelokopytova/202001051010polina_data/3DPredictor/input/H1/control.hic", "21", "21", "BP", 5000)
print(result)