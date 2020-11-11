from ChiPSeqReader import ChiPSeqReader

ctcf_reader = ChiPSeqReader("Z:/scratch/202002281332polina_data_2019/chip_mast_cells/mast_cells_WT/aquas_pipeline/peak/macs2/overlap/conservative_set/test3.bed")
ctcf_reader.read_file()
ctcf_reader.set_sites_orientation("Z:/scratch/202002281332polina_data_2019/3DPredictor/input/mast_cells/CTCF/SRR908255_mast_CTCF_orient.bed")
print(len(ctcf_reader.chr_data["chr7"].keys()))
print(ctcf_reader.chr_data["chr7"].keys())
ctcf_reader.chr_data["chr7"]["itemRGB"] = ["0,255,0"]*len(ctcf_reader.chr_data["chr7"])
print(ctcf_reader.chr_data["chr7"]["itemRGB"])
# ctcf_reader.chr_data["chr7"]["sigVal"] = ctcf_reader["plus_orientation"].apply(lambda x: x if )
plus_ori_ids = ctcf_reader.chr_data["chr7"][ctcf_reader.chr_data["chr7"]["plus_orientation"]!=0].index.tolist()
minus_ori_ids = ctcf_reader.chr_data["chr7"][ctcf_reader.chr_data["chr7"]["minus_orientation"]!=0].index.tolist()
print(minus_ori_ids)
print(ctcf_reader.chr_data["chr7"].loc[minus_ori_ids, 'itemRGB'])
ctcf_reader.chr_data["chr7"].loc[plus_ori_ids, 'itemRGB'] = "0,0,255"
ctcf_reader.chr_data["chr7"].loc[minus_ori_ids, 'itemRGB'] = "255,0,0"
print(ctcf_reader.chr_data["chr7"].loc[minus_ori_ids, 'itemRGB'])
# plus_ori_data = ctcf_reader.chr_data["chr7"].query("plus_orientation!=0")
# minus_ori_data = ctcf_reader.chr_data["chr7"].query("minus_orientation!=0")
# plus_ori_data[["chr", "start", "end"]].to_csv("D:/Users/Polina/3Dpredictor/input/K562/CTCF/chr7_plus_ori.bed", sep="\t", index=False)
# minus_ori_data[["chr", "start", "end"]].to_csv("D:/Users/Polina/3Dpredictor/input/K562/CTCF/chr7_minus_ori.bed", sep="\t", index=False)
ctcf_reader.chr_data["chr7"][["chr", "start", "end", "mids", "sigVal", "plus_orientation", "start", "end", "itemRGB"]].to_csv \
    ("Z:/scratch/202001051010polina_data/3DPredictor/input/H1/CTCF/H1_chr7_color.bed",sep="\t", index=False, header=False)
ctcf_reader.chr_data["chr7"][["chr","start","end","sigVal"]].to_csv \
    ("Z:/scratch/202001051010polina_data/3DPredictor/input/H1/CTCF/H1_chr7_color.bedgraph", sep="\t", index=False, header=False)