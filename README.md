--------
Contents
--------
This folder contains codes for implementing Half-approximate
matching. Please review the following paper for the serial
algorithm, based on Manne-Bisseling: 
http://www.staff.science.uu.nl/~bisse101/Articles/CP75.pdf

Also attached is a binary input file (karate.bin), for
Karate network 
( https://en.wikipedia.org/wiki/Zachary%27s_karate_club ).

This code requires an MPI library (MPI-3 compatible) 
and C++11 compliant compiler for building.

Please contact the following for any queries or support:

Sayan Ghosh, WSU (zsayanz at gmail dot com)
Mahantesh Halappanavar, PNNL (hala at pnnl dot gov)

-----
Cite
-----
Sayan Ghosh, Mahantesh Halappanavar, Ananth Kalyanaraman, Arif Khan, Assefaw Gebremedhin. 
Exploring MPI Communication Models for Graph Applications Using Graph Matching as a Case Study.
33rd IEEE International Parallel and Distributed Processing Symposium (IPDPS 2019).

-------
Compile
-------
Just invoking `make should build the program, without any
changes made to the Makefile. Please pass the appropriate
macro to select the MPI version, following are the macros
(also mentioned in the Makefile):

# Options for choosing the MPI variants
#-DUSE_MPI_NRM, RMA (passive mode) with neighbor communicator
#-DUSE_MPI_P2P, Plain P2P
#-DUSE_MPI_NCL, Neighborhood collective
#-DUSE_MPI_NPP, P2P with neighbor communicator
#-DUSE_MPI_RMA, Plain RMA (passive mode)
# Experimental:
#-DUSE_MPI_UPX, use -DREPLACE_UPX_WITH_RMA to replace UPCXX with MPI calls

-----
Input
-----
We require graphs to be converted from its native format to a binary format.
The binary converter is part of another application, please follow the 
instructions for using Vite for file conversion: https://github.com/Exa-Graph/vite

-------
Execute
-------
mpiexec -n 2 ./match -f karate.bin

Apart from using external file, it is possible to generate
in-memory a random geometric graph in a distributed fashion.

Possible options (can be combined):

1. -f <bin-file>   : Specify input binary file after this argument. 
2. -n <vertices>   : Pass total number of vertices of the generated graph.
3. -l              : Use distributed LCG for randomly choosing edges. If this option 
                     is not used, we will use C++ random number generator (using 
                     std::default_random_engine).
4. -p <percent>    : Specify percent of overall edges to be randomly generated between
                     processes.
5. -w              : Use Euclidean distance as edge weight. If this option is not used,
                     edge weights are considered as 1.0. Generate edge weight uniformly 
                     between (0,1) if Euclidean distance is not available (applicable to 
                     randomly generated edges).                    
6. -r <nranks>     : This is used to control the number of aggregators in MPI I/O and is
                     meaningful when an input binary graph file is passed with option "-f".
                     naggr := (nranks > 1) ? (nprocs/nranks) : nranks;

