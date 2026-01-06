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
#PBS -l select=1:ncpus=36:mem=120GB:mpiprocs=36
#PBS -l walltime=24:00:00
#PBS -P Personal
#PBS -j oe

cd $PBS_O_WORKDIR
NP=$(cat ${PBS_NODEFILE} | wc -l)
EOF

cat<<'EOF' >> $JobSub
##########    Load environment   ###########
module load intel/2022a
EOF

cat <<'EOF' >> $JobSub
##########    Activate conda     ###########
module load Miniconda3
conda activate ase
EOF

cat <<'EOF' >> $JobSub
##########    Set up VASP ENV    ###########
export VASP_COMMAND="mpirun -np ${NP} ${HOME}/opt/vasp/bin/vasp_632_vtst_solpp_beef_std"
export VASP_PP_PATH="${HOME}/opt/vasp/"
export ASE_VASP_VDW="${HOME}/opt/vasp/vdw/"

python run_vasp.py
EOF