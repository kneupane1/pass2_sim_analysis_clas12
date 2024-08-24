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
# Set HOME environment variable
export HOME=/home/kneupane

# Ensure TMP directory is writable
export TMPDIR=/scratch/slurm/$SLURM_JOB_ID/tmp
mkdir -p $TMPDIR

##export PATH=/work/clas12/users/tylern/software/hipo_tools/bin:$PATH
export PATH=/u/home/kneupane/hipo_tools_all_banks/hipo_tools/build/src/dst2root:$PATH
# Include path to parallel
#export PATH=/work/clas12/users/tylern/software:$PATH
export PATH=/u/home/kneupane/swif2_hipo2root:$PATH
# Verify environment setup
echo "======= Environment Setup ======="
echo "GCC_DIR: $GCC_DIR"
echo "ROOT_DIR: $ROOT_DIR"
echo "PATH: $PATH"
echo "LD_LIBRARY_PATH: $LD_LIBRARY_PATH"
which dst2root
which parallel
which hadd
echo "================================="

echo "======= dst2root ========="
#/work/clas12/users/tylern/software/parallel -j16 'dst2root -mc -a -t {} {.}.root' ::: *.hipo
/u/home/kneupane/swif2_hipo2root/parallel -j16 'dst2root -mc -t {} {.}.root' ::: *.hipo

if [ $? -ne 0 ]; then
    echo "Error: dst2root failed"
    exit 1
fi

# Verify .root files creation
num_root_files=$(ls -1 *.root 2>/dev/null | wc -l)
if [ "$num_root_files" -eq 0 ]; then
    echo "Error: No .root files were created"
    exit 1
fi
echo "Created $num_root_files .root files"
echo "==========================="

echo "======= hadd ========="
# Create temporary directory
mkdir -p tmp

# Merge .root files into merged.root
##/cvmfs/oasis.opensciencegrid.org/jlab/scicomp/sw/el9/root/6.30.04-gcc11.4.0/bin/hadd -j 16 -d $PWD/tmp merged.root *.root
/cvmfs/oasis.opensciencegrid.org/jlab/scicomp/sw/el9/root/6.30.04-gcc11.4.0/bin/hadd -k -j 16 -d $PWD/tmp merged.root *.root

if [ $? -ne 0 ]; then
    echo "Error: hadd failed"
    exit 1
fi


# Check if merged.root was created
if [ ! -f merged.root ]; then
    echo "Error: merged.root not found"
    exit 1
fi
echo "merged.root created successfully"
echo "==========================="

