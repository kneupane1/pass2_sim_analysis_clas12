#!/bin/bash


# Manually set up GCC environment
export GCC_DIR=/cvmfs/oasis.opensciencegrid.org/jlab/scicomp/sw/el9/gcc/11.4.0
export PATH=$GCC_DIR/bin:$PATH
export LD_LIBRARY_PATH=$GCC_DIR/lib64:$LD_LIBRARY_PATH

# Print the GCC version to confirm the correct setup
gcc --version

# Manually set up ROOT environment (update the path as needed)
export ROOT_DIR=/cvmfs/oasis.opensciencegrid.org/jlab/scicomp/sw/el9/root/6.30.04-gcc11.4.0
export PATH=$ROOT_DIR/bin:$PATH
export LD_LIBRARY_PATH=$ROOT_DIR/lib:$LD_LIBRARY_PATH

# Current date
_now=$(date +"%m_%d_%Y")

##source /group/clas12/packages/setup.sh

##module load gcc/9.2.0
##module load root/6.26.10

export PATH=/u/home/kneupane/sim_multi_thr_clas12_ana/build:$PATH


echo "======= clas12_mc ========="
clas12_mc out.root *.root
echo "======= clas12_mc ========="

