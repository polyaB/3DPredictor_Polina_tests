#PBS -q vkop2q
#PBS -l select=1:ncpus=1:mem=150g
#PBS -l walltime=24:00:00
#PBS -l cput=240:00:00

cd ~/3DPredictor_Polina_tests/; module load jre/1.8.0; python train_and_validate_K562_feature_imp_test.py