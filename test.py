
import pandas as pd
import numpy as np
data1 = pd.read_csv("Z:/scratch/202011291709Polya_data/3DPredictor/input/mast_cells/RNA-seq/Our_rna_seq/wt_1.fq.tabular.pre.sorted.bedgraph", sep="\t",
                   names=["chr", "start", "end", "FPKM"])
data2 = pd.read_csv("Z:/scratch/202011291709Polya_data/3DPredictor/input/mast_cells/RNA-seq/Our_rna_seq/wt_2.fq.tabular.pre.sorted.bedgraph", sep="\t",
                   names=["chr", "start", "end", "FPKM"])
data3 = pd.read_csv("Z:/scratch/202011291709Polya_data/3DPredictor/input/mast_cells/RNA-seq/Our_rna_seq/wt_3.fq.tabular.pre.sorted.bedgraph", sep="\t",
                   names=["chr", "start", "end", "FPKM"])
r = np.corrcoef(data1["FPKM"], data2["FPKM"])
r2 = np.corrcoef(data1["FPKM"], data3["FPKM"])
print(r)
print(r2)
# data["chr"] = data["chr"].apply(lambda x: "chr"+str(x))
# data.to_csv("Z:/scratch/202011291709Polya_data/3DPredictor/input/mast_cells/RNA-seq/Our_rna_seq/wt_1.fq.tabular.pre.sorted2.bedgraph",
#             sep="\t", index=False, header=False)
