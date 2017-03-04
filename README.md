# Parallel_CAVIAR
In current version of CAVIAR (https://github.com/fhormoz/caviar), it runs in a serial manner. To speed up the running, I convert it to a parallel version, which applies openmp to run with multi-thread. Simulation results demonstrate that when explore variants are massive, the running time is only 1/nthread of original version. 
The main changes happen in PostCal.cpp. In previous version, when calculating the likelihood of the causal states, it was run one by one. In this parallel version, firstly, I store the causal status in a 2-dimensional vector, and then visit causal status in parallel. 
The parallel version can be much faster when multiple CPUs are available. To run it, please add "-t n", n is the number of threads. 
