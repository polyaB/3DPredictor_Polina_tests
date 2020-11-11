from matplot2hic import MatPlot2HiC
from matrix_plotter import MatrixPlotter
import pandas as pd
# "/mnt/scratch/ws/psbelokopytova/201903101031polina/3DPredictor/out/analysis/rearrangements/del_unc5b/WT_hepat_predicted_new.txt"
# "/mnt/scratch/ws/psbelokopytova/201903101031polina/3DPredictor/out/analysis/rearrangements/del_unc5b/del_hepat_predicted_final.txt"
# data = pd.read_csv("/mnt/scratch/ws/psbelokopytova/201903101031polina/3DPredictor/out/analysis/rearrangements/del_unc5b/del_thymus_predicted.txt", sep="\t")
# print(data)
# data["chr"] = data["chr"].apply(lambda  x: str(x))
#
# control_data = pd.read_csv("/mnt/scratch/ws/psbelokopytova/201903101031polina/3DPredictor/out/analysis/rearrangements/del_unc5b/WT_hepat.txt", sep="\t")
# control_data["chr"] = control_data["chr"].apply(lambda  x: str(x))
# mp=MatrixPlotter()
# mp.set_data(data)
# mp.set_control(control_data)
# MatPlot2HiC(mp,"thymus_predDel-WTHepat_unc5b", "/mnt/scratch/ws/psbelokopytova/201903101031polina/3DPredictor/out/analysis/rearrangements/")

# data_Rowley = pd.read_csv("/mnt/scratch/ws/psbelokopytova/201911101015polinas_data/3DPredictor/out/scc/chr4-model-both-ctcf-and-gro.norm_contacts.scc", sep=" ")
# data_Rowley["chr"]=["chr4"]*len(data_Rowley)
# data_Pierro = pd.read_csv("/mnt/scratch/ws/psbelokopytova/201911101015polinas_data/3DPredictor/input/GM12878/chr4.GM12878.DiPierro.50KB.contacts",
#                           sep="\t", names=["contact_st", "contact_en", "contact_count"])
data_Qi = pd.read_csv("/mnt/scratch/ws/psbelokopytova/201911101015polinas_data/3DPredictor/out/scc/forscc_sparse_merge3DPred_contact_map_comb_chr1_5kb.txt",
                          sep=" ")
data_Qi["chr"]=["chr1"]*len(data_Qi)
# data_Pierro["chr"] = ["chr4"]*len(data_Pierro)
# data=data_Pierro[["chr","contact_st", "contact_en", "contact_count"]]
# control_data=data_Pierro[["chr", "contact_st", "contact_en", "contact_count"]]
# print(data_Rowley.keys())
# data = data_Rowley[["chr", "contact_st", "contact_en", "0"]]
# data.rename(columns={"0":"contact_count"}, inplace=True)
# control_data=data_Rowley[["chr", "contact_st", "contact_en", "0"]]
# control_data.rename(columns={"0":"contact_count"}, inplace=True)
# data = pd.read_csv("/mnt/scratch/ws/psbelokopytova/201911101015polinas_data/3DPredictor/out/scc/add_random_lig.chr4.5000.model8669235584.validation..equal.GM12878.1124356093.scc",
#                      sep="\t")
# data["chr"]=["chr4"]*len(data)
mp=MatrixPlotter()
# mp.set_data(data)
# mp.set_control(control_data)
# mp.set_control(control_data[["chr", "contact_st", "contact_en", "contact_count"]])
# print(data[["chr", "contact_st", "contact_en", "contact_count"]])
# print(data[["chr", "contact_st", "contact_en", "0"]])
# MatPlot2HiC(mp,"Rowley_pred_WT", "/mnt/scratch/ws/psbelokopytova/201911101015polinas_data/3DPredictor/out/hic_files/")

data_3Dpredictor = pd.read_csv("/mnt/scratch/ws/psbelokopytova/201911101015polinas_data/3DPredictor/out/scc/add_random_chr1.5000.model8669235584.validation..equal.GM12878.3108998743.scc",
                     sep=" ")
print(data_3Dpredictor.keys())
# data_experiment =pd.read_csv("/mnt/scratch/ws/psbelokopytova/201911101015polinas_data/3DPredictor/input/GM12878/chr4.50KB.GM12878.control.contacts",
#                                sep="\t", names=["contact_st", "contact_en", "contact_count"])
# data_3Dpredictor = pd.read_csv("/mnt/scratch/ws/psbelokopytova/201911101015polinas_data/3DPredictor/input/GM12878/chr4.50KB.GM12878.3dpredictor.contacts",
#                                sep="\t", names=["contact_st", "contact_en", "contact_count"])
# data_experiment["chr"]=["chr4"]*len(data_experiment)
data_3Dpredictor["chr"]=["chr1"]*len(data_3Dpredictor)
# print(data_3Dpredictor.keys())
data = data_3Dpredictor[["chr", "contact_st", "contact_en", "0"]]
data.rename(columns={"0":"contact_count"}, inplace=True)
control_data=data_Qi[["chr", "contact_st", "contact_en", "contact_count"]]
# control_data.rename(columns={"0":"contact_count"}, inplace=True)
# data["contact_count"] = data["contact_count"]*11478
# pierro_3D = pd.merge(data_Pierro,data_experiment, on=["chr", "contact_st", "contact_en"], how="inner")
# pierro_3D = pd.merge(data_3Dpredictor,data_experiment, on=["chr", "contact_st", "contact_en"], how="inner")
#
# print("232323232323232323232322323232323323232323232keys", pierro_3D.keys())
# data=pierro_3D[["chr", "contact_st", "contact_en", "contact_count_x"]]
# pierro_3D.rename(columns={"contact_count_y":"contact_count", "contact_count_x":"0"}, inplace=True)
# control_data = pierro_3D[["chr", "contact_st", "contact_en", "contact_count_y"]]
# control_data.rename(columns={"contact_count_y":"contact_count"}, inplace=True)
# pierro_3D.rename(columns={"contact_count_y":"contact_count", "contact_count_x":"0"}, inplace=True)
# pierro_3D[["contact_st", "contact_en", "contact_count", "0"]].to_csv("/mnt/scratch/ws/psbelokopytova/201911101015polinas_data/3DPredictor/out/scc/experiment_3D_50KB.scc", sep=" ", index=False)
# print("done")
# control_data.rename(columns={"contact_count_y":"contact_count"}, inplace=True)
# print(min(control_data["contact_count"]), max(control_data["contact_count"]), min(data["contact_count"]), max(data["contact_count"]))
# control_data=data_experiment[["chr", "contact_st", "contact_en", "contact_count"]]
mp.set_control(control_data[["chr", "contact_st", "contact_en", "contact_count"]])
mp.set_data(data[["chr", "contact_st", "contact_en", "contact_count"]])
MatPlot2HiC(mp,"add_random_3DPredictor_experiment_5KB", "/mnt/scratch/ws/psbelokopytova/201911101015polinas_data/3DPredictor/out/hic_files/")
# MatPlot2HiC(mp,"add_random_3DPredictor_Qi_5KB", "/mnt/scratch/ws/psbelokopytova/201911101015polinas_data/3DPredictor/out/hic_files/")
