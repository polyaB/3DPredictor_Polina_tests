import pandas as pd
from matrix_plotter import MatrixPlotter
from matplot2hic import MatPlot2HiC

out_dir = "/mnt/scratch/ws/psbelokopytova/201901151331psbelokopytova/3DPredictor/out/hic_files"
fname = "Hepat-Hepat_Hepat-NPC"
predicted_data = pd.read_csv("/mnt/scratch/ws/psbelokopytova/201901151331psbelokopytova/3DPredictor/out/pics/scc/hepat-NPC/Hepatpr_N.txt", sep="\t")
control_data= pd.read_csv("/mnt/scratch/ws/psbelokopytova/201901151331psbelokopytova/3DPredictor/out/pics/scc/hepat-NPC/Hepatpr_h.txt", sep="\t")
mp = MatrixPlotter()
mp.set_data(predicted_data)
mp.set_control(control_data)
MatPlot2HiC(mp,fname, out_dir)
