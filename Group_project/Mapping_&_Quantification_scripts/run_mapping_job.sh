env2lmod
module load minimap2/2.17
module load bedtools2/2.30.0
conda activate sta
sbatch --cpus-per-task=8 --mem-per-cpu=4G --wrap="source $HOME/STA/run_mapping_05.sh" -o $HOME/STA/job_05.out
sbatch --cpus-per-task=8 --mem-per-cpu=4G --wrap="source $HOME/STA/run_mapping_06.sh" -o $HOME/STA/job_06.out --time 3:00:00
sbatch --cpus-per-task=8 --mem-per-cpu=4G --wrap="source $HOME/STA/run_mapping_07.sh" -o $HOME/STA/job_07.out --time 3:00:
sbatch --cpus-per-task=8 --mem-per-cpu=4G --wrap="source $HOME/STA/run_mapping_08.sh" -o $HOME/STA/job_08.out