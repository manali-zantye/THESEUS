#!/bin/bash
##ENVIRONMENT SETTINGS; CHANGE WITH CAUTION
#SBATCH --export=NONE        #Do not propagate environment
#SBATCH --get-user-env=L     #Replicate login environment
  
##NECESSARY JOB SPECIFICATIONS
#SBATCH --job-name=thses     #Set the job name
#SBATCH --time=30:00:00            #Set the wall clock limit
#SBATCH --ntasks=1                 #Request 1 task
#SBATCH --ntasks-per-node=1        #Request 1 task/core per node
#SBATCH --mem=10000M                #Request 2560MB (2.5GB) per node
#SBATCH --output=Output_file.%j.out    #Send stdout/err to "Example1Out.[jobID]"

#SBATCH --mail-type=ALL              #Send email on all job events
#SBATCH --mail-user=manali.137@tamu.edu    #Send all emails to email_address 

ml Anaconda3
ml gams



gams battery_test.gms u1 = 145 u2 = 0.1667 u3 = 24