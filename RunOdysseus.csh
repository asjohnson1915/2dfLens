#PBS -S /bin/csh
#PBS -l walltime=10:00:00
#PBS -l nodes=8:ppn=1
#PBS -l pmem=10gb
#PBS -q sstar
#PBS -d /lustre/projects/p018_swin/ajohnson/2dfLens
#PBS -N calc_shear_galaxy
#PBS -o logs/2df_shear.log
#PBS -e logs/2df_shear.err 

module load python/2.7.2
module load intel/13.3.1
module unload openmpi
module load mpi4py

mpiexec -n 8 python Odysseus.py > logs/2df_shear.txt