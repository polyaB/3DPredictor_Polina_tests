#PBS -q xl230g9q
#PBS -l select=1:ncpus=1:mem=50g
#PBS -l walltime=3:00:00
#PBS -l cput=30:00:00

cd ~/3DPredictor_Polina_tests/; module load jre/1.8.0; python train_and_validate_K562.py