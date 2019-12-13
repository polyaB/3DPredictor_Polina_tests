#PBS -q xl230g9q 
#PBS -l select=1:ncpus=12:mem=60g
#PBS -l walltime=48:00:00
#PBS -l cput=480:00:00
#PBS -k oe

cd ~/3DPredictor_Polina_tests/;python generata_data_K562_cluster.py 14 'contacts.gz'

