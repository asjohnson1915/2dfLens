#PBS -S /bin/csh
#PBS -l walltime=01:00:00
#PBS -l nodes=1:ppn=1
#PBS -l mem=1gb
#PBS -q sstar
#PBS -d /lustre/projects/p018_swin/ajohnson/2dfLens
#PBS -N 2df_SortSims
#PBS -o logs/SortSimsNew.log
#PBS -e logs/SortSimsNew.err 

module load python/2.7.2
module unload openmpi
module load mpi4py

mpiexec -n 1 python SortMocks.py > logs/SortSimsNew.txt