# pass2_analysis_clas12

This is for the two-pion channel analysis using clas12 rg-a data converted to the root files using dst2root
Prerequisites:

- [cern root](https://root.cern.ch/)

### cpp

A c++ example can be found in the cpp folder.

To build:

```bash
make clean && make
```

To run:

```bash
CLAS12_E=10.6041 NUM_THREADS=4 ./clas12_analysis output.root /path/to/input/*.root
```

If it breaks, reduce the number of threads for the number of files. In general each thread should have 2 or more files and the number of threads should be less than or equal to the number of cores you are using.

So for 16 files on a 4 core computer use NUM_THREADS=4 and each thread will process 4 files. For 4 files on a 4 core computer use NUM_THREADS=1 or NUM_THREADS=2 so each thread will have 4 or 2 files respecfully.
