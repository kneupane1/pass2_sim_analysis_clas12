#!/bin/bash

# Loop over Q² values from 1 to 8
for ((q2=1; q2<=8; q2++)); do
    echo "Compiling mPim_cs_50_MeV_W_bin.cxx..."
    g++ mPim_cs_50_MeV_W_bin.cxx -o cs_measure `root-config --cflags --libs`
    
    if [ $? -ne 0 ]; then
        echo "Compilation failed. Exiting..."
        exit 1
    fi

    echo "Running cs_measure for Q² = $q2..."
    ./cs_measure $q2  # Pass Q² as an argument to the program
    
    echo "Finished Q² = $q2"
    echo "--------------------------------"
done

echo "All Q² values processed!"

