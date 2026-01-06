#!/bin/bash

# This script prepares the environment and submit script for vasp calculations
# on the NSCC cluster. It sets up necessary environment variables and creates
# a job submission script compatible with the cluster's job scheduler.

Nprocs=16
RootDir=`pwd |rev|awk -F "/" '{print $1}'|rev`

# Structure
VaspStru=$1
WorkDir=`echo $VaspStru | sed 's/.vasp//g'`
WorkerID=`echo $WorkDir | sed 's/.*_\([0-9]\+\)$/\1/'`
JobSub="sub_${WorkerID}.sh"

if [ ! -d $WorkDir ]; then
    mkdir $WorkDir
fi
cd $WorkDir
mv ../$VaspStru POSCAR

cat<<'EOF' > $JobSub
#!/bin/sh -f
#PBS -N vasp_job_${WorkerID}
#PBS -l select=1:ncpus=64:mem=120GB:mpiprocs=64
#PBS -l walltime=24:00:00
#PBS -q normal
#PBS -P Personal
#PBS -j oe

cd $PBS_O_WORKDIR
NP=$(cat ${PBS_NODEFILE} | wc -l)
EOF

cat<<'EOF' >> $JobSub
##########    Load environment   ###########
module swap craype-x86-rome craype-x86-milan
module swap PrgEnv-cray PrgEnv-gnu
module load cray-hdf5-parallel/1.12.1.1
module swap craype-network-ofi craype-network-ucx
module swap cray-mpich/8.1.15 cray-mpich-ucx/8.1.15
module load intel
module load gcc
module load mkl
EOF

cat <<'EOF' >> $JobSub
##########    Activate conda     ###########
module load miniforge3/23.10
conda activate ase
export PYTHONPATH=/home/users/nus/c_yang/opt/ase:$PYTHONPATH
EOF

cat <<'EOF' >> $JobSub
##########    Set up VASP ENV    ###########
export VASP_COMMAND="mpirun -np ${NP} ${HOME}/project-1/opt/vasp/bin/vasp_651_intel_mkl_std"
export VASP_PP_PATH="${HOME}/opt/vasp/pseudopotetnial"
export ASE_VASP_VDW="${HOME}/opt/vasp/vdw/"

python run_vasp.py
EOF