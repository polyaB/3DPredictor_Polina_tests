__version__ = "0.2beta"

import logging
import pandas as pd
import sys
import os

sourcedir = os.path.join(os.path.dirname(os.path.abspath(sys.argv[0])), "../3Dpredictor/source")
sys.path.append(sourcedir)
from ChiPSeqReader import ChiPSeqReader
import pickle
from PredictorGenerators import SmallChipSeqPredictorGenerator, SitesOrientPredictorGenerator, \
    OrientBlocksPredictorGenerator, ConvergentPairPredictorGenerator
from shared import Interval, Parameters
from RNASeqReader import RNAseqReader
from DataGenerator import DataGenerator
from check_file_formats import *
import argparse
import datetime

logging.basicConfig(format='%(asctime)s %(message)s', datefmt='%I:%M:%S', level=logging.DEBUG)

if __name__ == '__main__':
    RNA_seq_file = "/mnt/scratch/ws/psbelokopytova/202011291709Polya_data/3DPredictor/input/mast_cells/RNA-seq/Our_rna_seq/StringTie on collection 64_ Gene abundance estimates/wt_1.fq.tabular"
    CTCF_file = "/mnt/scratch/ws/psbelokopytova/202011291709Polya_data/3DPredictor/input/NPC/CTCF/GSE96107_NPC_CTCF.IDR0.05.filt.narrowPeak"
    CTCF_orient_file = "/mnt/scratch/ws/psbelokopytova/202011291709Polya_data/3DPredictor/input/NPC/CTCF/GSE96107_NPC_CTCF.IDR0.05.filt.narrowPeak-orient.bed"
    chr = sys.argv[1]
    conttype = "contacts.gz"
    interval_start = sys.argv[2]
    interval_end = sys.argv[3]
    resolution = 5000
    # model_path = "/mnt/scratch/ws/psbelokopytova/202011291709Polya_data/3DPredictor/out/models/912470444"
    model_path = "/mnt/scratch/ws/psbelokopytova/202011291709Polya_data/3DPredictor/out/models/models_for_web_3DPredictor/mouse_NPC_5KB_oddChrms_CtcfRnaDist"
    out_file = "/mnt/scratch/ws/psbelokopytova/202011291709Polya_data/3DPredictor/out/mast_cells/"+chr+"_"+str(interval_start)+"_"+str(interval_end)+\
               "_"+str(resolution)+model_path.split("/")[-1]+"_prediction"

    params = Parameters()
    params.window_size = int(resolution)  # region around contact to be binned for predictors
    params.mindist = int(resolution) * 2 + 1  # minimum distance between contacting regions
    params.maxdist = 1500000  # maximum distance between contacting regions
    params.max_cpus = 11
    # params.keep_only_orient=False #set True if you want use only CTCF with orient
    params.multiprocessing = True
    # params.write_to_file = False
    # Read CTCF data
    # CTCF_file format: ENCODE narrow peak
    # CTCF_orient_file format: chr--start--end--name--score--strand
    logging.info('create CTCF_PG')
    # set path to the CTCF chip-seq file:
    params.ctcf_reader = ChiPSeqReader(CTCF_file, name="CTCF")
    params.ctcf_reader.read_file()
    # set path to the CTCF_orient file:
    params.ctcf_reader.set_sites_orientation(CTCF_orient_file)
    # set corresponding predictor generators and its options:
    OrientCtcfpg = SitesOrientPredictorGenerator(params.ctcf_reader, N_closest=4)
    NotOrientCTCFpg = SmallChipSeqPredictorGenerator(params.ctcf_reader, params.window_size, N_closest=4)
    ctcf_reader_orientOnly = ChiPSeqReader(CTCF_file, name="CTCF")
    ctcf_reader_orientOnly.read_file()
    ctcf_reader_orientOnly.set_sites_orientation(CTCF_orient_file)
    ctcf_reader_orientOnly.keep_only_with_orient_data()
    # set corresponding predictor generators and its options:
    OrientBlocksCTCFpg = OrientBlocksPredictorGenerator(ctcf_reader_orientOnly, params.window_size)
    ConvergentPairPG = ConvergentPairPredictorGenerator(params.ctcf_reader, binsize=params.window_size)

    # Read RNA-Seq data
    # RNA-seq_file format: this file should have fields "gene", "start", "end", "chr","FPKM"
    # you can rename table fields below
    params.RNAseqReader = RNAseqReader(RNA_seq_file, name="RNA")
    # read RNA-seq data and rename table fields
    params.RNAseqReader.read_file(rename={ "Gene ID": "gene",
                              "Start": "start",
                              "End": "end",
                              "Reference": "chr",
                              "FPKM": "sigVal"}, sep="\t")
    # set corresponding predictor generators and its options:
    RNAseqPG = SmallChipSeqPredictorGenerator(params.RNAseqReader, window_size=params.window_size, N_closest=3)
    params.pgs =  [OrientCtcfpg, NotOrientCTCFpg, OrientBlocksCTCFpg,ConvergentPairPG,RNAseqPG]

    interval_start_bin = int(interval_start) // int(resolution) * int(resolution)
    interval_end_bin = int(interval_end) // int(resolution) * int(resolution)

    logging.info("dump trained model")
    trained_predictor = pickle.load(open(model_path, "rb"))

    outfile = open(out_file, "w")
    outfile.write("chr" + "\t" + "contact_st" + "\t" + "contact_en" + "\t" + "predicted_contact_count" + "\n")

    for contact_st in range(interval_start_bin, interval_end_bin + resolution, int(resolution)):
        logging.info(str(datetime.datetime.now()) + " " + str(contact_st))
        for contact_en in range(contact_st, interval_end_bin + resolution, int(resolution)):
            if (contact_en - contact_st) <= params.maxdist and (contact_en - contact_st) >= params.mindist:
                contact_df = pd.DataFrame([[chr, contact_st, contact_en, contact_en - contact_st]],
                                          columns=['chr', 'contact_st', 'contact_en', 'dist'])
                generator = DataGenerator()
                result = generator.contact2predictors(contact_df, params)
                result_df = pd.DataFrame([[chr, contact_st, contact_en, contact_en - contact_st, 0] + result[1]],
                                         columns=['chr', 'contact_st', 'contact_en', 'contact_dist', "contact_count"] +
                                                 result[0])
                trained_predictor.validate(result_df, show_plot=False, out_dir=out_file, df_input=True)
                # print("predicted")
                # print(trained_predictor.predicted)
                outfile.write(str(chr) + "\t" + str(contact_st) + "\t" + str(contact_en) + "\t" + str(
                    trained_predictor.predicted[0]) + "\n")
                #     break
        # break
    outfile.close()
