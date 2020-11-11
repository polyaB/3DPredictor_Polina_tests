import pandas as pd

data_GM12878 = pd.read_csv("Z:/scratch/202010021628data_Polina/3DPredictor/out/scc/Qi_for_scc_all_contacts.scc",
                   sep=" ")
# data_K562 = pd.read_csv("Z:/scratch/202010021628data_Polina/3DPredictor/out/scc/chr1.5000.model3572004608.validation..equal.K562.3108998743.scc",
#                         sep=" ")
# data_Gm12878_rep1 = pd.read_csv("D:/Users/Polina/my_papers/3D_review/figures_for_review/GM12878_007_rep_contacts.txt",
#                                 sep="\t", names=["contact_st", "contact_en", "count"])
# print(data_Gm12878_rep1.keys())
# data_Gm12878_rep2 = pd.read_csv("D:/Users/Polina/my_papers/3D_review/figures_for_review/GM12878_003_rep_contacts.txt",
#                                 sep="\t", names=["contact_st", "contact_en", "count"])
data_IMR = pd.read_csv("D:/Users/Polina/my_papers/3D_review/IMR90_combined_contacts.txt", sep="\t", names=["contact_st", "contact_en", "count"])



# print(data_GM12878.keys())
# merge_data = data_GM12878.merge(data_K562, how='left', on=["contact_st", "contact_en"], indicator=True)
# merge_data = data_Gm12878_rep1.merge(data_Gm12878_rep2, how='left', on=["contact_st", "contact_en"])
merge_data_2 = data_GM12878.merge(data_IMR, how='left', on=["contact_st", "contact_en"], indicator=True)
print(merge_data_2.keys())

new_data = merge_data_2[merge_data_2["_merge"]=="both"]
# new_data[["contact_st", "contact_en", "count_x", "count_y"]].to_csv("D:/Users/Polina/my_papers/3D_review/figures_for_review/GM12878_K562.scc",
#                                                           sep=" ", header=["contact_st", "contact_en", "contact_count", "0"])
new_data.dropna(inplace=True)
new_data[["contact_st", "contact_en", "count", "0"]].to_csv("D:/Users/Polina/my_papers/3D_review/figures_for_review/IMR90.scc",
                                                          sep=" ", header=["contact_st", "contact_en", "contact_count", "0"])
# print(merge_data)