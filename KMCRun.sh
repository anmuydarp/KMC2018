#!/bin/bash  

#SBATCH -n 1 -t 1-12:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=pmuralid@asu.edu

module load gcc/7.3.0
                                                                                                                                                                       
/home/pmuralid/KMC2018/Code/KineticMC
