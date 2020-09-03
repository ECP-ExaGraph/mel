#include <sys/resource.h>
#include <sys/time.h>
#include <unistd.h>

#include <cassert>
#include <cstdlib>

#include <fstream>
#include <iostream>
#include <sstream>
#include <string>

#include <mpi.h>

#if defined(USE_MPI_RMA)
#include "maxematch_rma.hpp"
#elif defined(USE_MPI_NCL)
#include "maxematch_ncl.hpp"
#elif defined(USE_MPI_NCN)
#include "maxematch_ncn.hpp"
#elif defined(USE_MPI_NRM)
#include "maxematch_nrm.hpp"
#elif defined(USE_MPI_NPP)
#include "maxematch_npp.hpp"
#elif defined(USE_MPI_P2P)
#include "maxematch_p2p.hpp"
#elif defined(USE_MPI_UPX)
#include "maxematch_upx.hpp"
#else
#include "maxematch_serial.hpp"
#endif

static std::string inputFileName;
static int me, nprocs;
static int ranksPerNode = 1;
static GraphElem nvRGG = 0;
static int generateGraph = 0;
static int randomEdgePercent = 0;
static bool randomNumberLCG = false;
static bool readBalanced = false;
static double threshold = 1.0E-6;
static bool isUnitEdgeWeight = true;

// parse command line parameters
static void parseCommandLine(const int argc, char * const argv[]);

int main(int argc, char *argv[])
{
    double t0, t1, td, td0, td1;
    int max_threads;
    
    max_threads = omp_get_max_threads();
    
    if (max_threads > 1) {
      int provided;
      MPI_Init_thread(&argc, &argv, MPI_THREAD_FUNNELED, &provided);
      if (provided < MPI_THREAD_FUNNELED) {
          std::cerr << "MPI library does not support MPI_THREAD_FUNNELED." << std::endl;
          MPI_Abort(MPI_COMM_WORLD, -99);
      }
    } else {
	    MPI_Init(&argc, &argv);
    }

    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &me);
  
#if defined(USE_MPI_UPX)
    upcxx::init();
#endif

    // command line options
    parseCommandLine(argc, argv);
 
    MPIGraph* g = nullptr;
    
    td0 = MPI_Wtime();

    // generate graph only supports RGG as of now
    if (generateGraph) 
    { 
        GenerateRGG gr(nvRGG);
        g = gr.generate(randomNumberLCG, isUnitEdgeWeight, randomEdgePercent);

        //g->print(false);

        if (me == 0) 
        {
            std::cout << "Generated Random Geometric Graph with d: " << gr.get_d() << std::endl;
            const GraphElem nv = g->get_nv();
            const GraphElem ne = g->get_ne();
            std::cout << "Number of vertices: " << nv << std::endl;
            std::cout << "Number of edges: " << ne << std::endl;
            //std::cout << "Sparsity: "<< (double)((double)nv / (double)(nvRGG*nvRGG))*100.0 <<"%"<< std::endl;
            std::cout << "Average degree: " << (ne / nv) << std::endl;
        }

        MPI_Barrier(MPI_COMM_WORLD);
    }
    else // read input graph
    {
        BinaryEdgeList rm;
        if (readBalanced == true)
            g = rm.read_balanced(me, nprocs, ranksPerNode, inputFileName);
        else
            g = rm.read(me, nprocs, ranksPerNode, inputFileName);
        //g->print();
    }
        
    g->print_dist_stats();
    assert(g != nullptr);

    MPI_Barrier(MPI_COMM_WORLD);
#ifdef DEBUG_PRINTF  
    assert(g);
#endif
    td1 = MPI_Wtime();
    td = td1 - td0;

    double tdt = 0.0;
    MPI_Reduce(&td, &tdt, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

    if (me == 0)  
    {
        if (!generateGraph)
            std::cout << "Time to read input file and create distributed graph (in s): " 
                << tdt << std::endl;
        else
            std::cout << "Time to generate distributed graph of " 
                << nvRGG << " vertices (in s): " << tdt << std::endl;
    }

    // start graph matching

#if defined(USE_MPI_RMA)
    if (me == 0)
        std::cout << "MPI-3 RMA: ";
    MaxEdgeMatchRMA mt(g); 
    mt.create_mpi_win(); // create MPI windows
#elif defined(USE_MPI_NCL)
    if (me == 0)
        std::cout << "MPI-3 Neighborhood Collectives: ";   
    MaxEdgeMatchNCL mt(g);
#elif defined(USE_MPI_NCN)
    if (me == 0)
        std::cout << "MPI-3 Neighborhood Nonblocking Collectives: ";   
    MaxEdgeMatchNCN mt(g);
#elif defined(USE_MPI_NRM)
    if (me == 0)
        std::cout << "MPI-3 Neighborhood RMA: ";   
    MaxEdgeMatchNRM mt(g);
#elif defined(USE_MPI_NPP)
    if (me == 0)
        std::cout << "MPI-3 Neighborhood P2P: ";   
    MaxEdgeMatchNPP mt(g);
#elif defined(USE_MPI_P2P)
    if (me == 0)
        std::cout << "MPI-2 Nonblocking Send/Recv: ";
    MaxEdgeMatchP2P mt(g);
#elif defined(USE_MPI_UPX)
    if (me == 0)
        std::cout << "MPI and UPCXX: ";
    MaxEdgeMatchUPX mt(g);      
#if defined(REPLACE_UPX_WITH_RMA)
    mt.create_mpi_win(); // create MPI windows
#else
    // create dist_object<...> per process
    upcxx::dist_object<upcxx::global_ptr<GraphElem>> dobj(upcxx::new_array<GraphElem>(mt.get_nelems()));
    mt.create_upx_ptr(dobj);
#endif
#else
    std::cout << "Serial: ";
    MaxEdgeMatch mt(g);

    if (nprocs > 1) // error
    {
        if (me == 0)
            std::cout << "Multiple processes detected, choose a parallel matching implementation. Exiting..." << std::endl;
        MPI_Abort(MPI_COMM_WORLD, -99);
    }
#endif

    MPI_Barrier(MPI_COMM_WORLD);
    
    // invoke matching
    t0 = MPI_Wtime();
    std::vector<EdgeTuple> const& M = mt();
    
    MPI_Barrier(MPI_COMM_WORLD);
    
    t1 = MPI_Wtime();
    double p_tot = t1 - t0, t_tot = 0.0;
    
    MPI_Reduce(&p_tot, &t_tot, 1, MPI_DOUBLE, 
            MPI_SUM, 0, MPI_COMM_WORLD);
    if (me == 0)
        std::cout << "Execution time (in s) for maximal edge matching: " 
            << (double)(t_tot/(double)nprocs) << std::endl;

#if defined(CHECK_RESULTS)    
    mt.check_results();
#endif
#if defined(PRINT_RESULTS)    
    mt.print_M();
#endif
 
    MPI_Barrier(MPI_COMM_WORLD);
    
#if defined(USE_MPI_UPX)
#if defined(REPLACE_UPX_WITH_RMA)
    mt.destroy_mpi_win();
#else
    mt.destroy_upx_ptr();
    upcxx::finalize();
#endif
#endif

    // destroy win
#if defined(USE_MPI_RMA)
    mt.destroy_mpi_win();
#endif
    mt.clear();

    MPI_Finalize();

    return 0;
}

void parseCommandLine(const int argc, char * const argv[])
{
  int ret;

  while ((ret = getopt(argc, argv, "f:br:n:wlp:")) != -1) {
    switch (ret) {
    case 'f':
      inputFileName.assign(optarg);
      break;
    case 'b':
      readBalanced = true;
      break;
    case 'r':
      ranksPerNode = atoi(optarg);
      break;
    case 'n':
      nvRGG = atol(optarg);
      if (nvRGG > 0)
          generateGraph = true; 
      break;
    case 'w':
      isUnitEdgeWeight = false;
      break;
    case 'l':
      randomNumberLCG = true;
      break;
    case 'p':
      randomEdgePercent = atoi(optarg);
      break;
    default:
      assert(0 && "Should not reach here!!");
      break;
    }
  }

  if (me == 0 && (argc == 1)) {
      std::cerr << "Must specify some options." << std::endl;
      MPI_Abort(MPI_COMM_WORLD, -99);
  }
  
  if (me == 0 && !generateGraph && inputFileName.empty()) {
      std::cerr << "Must specify a binary file name with -f or provide parameters for generating a graph." << std::endl;
      MPI_Abort(MPI_COMM_WORLD, -99);
  }
   
  if (me == 0 && !generateGraph && randomNumberLCG) {
      std::cerr << "Must specify -g for graph generation using LCG." << std::endl;
      MPI_Abort(MPI_COMM_WORLD, -99);
  } 
   
  if (me == 0 && !generateGraph && randomEdgePercent) {
      std::cerr << "Must specify -g for graph generation first to add random edges to it." << std::endl;
      MPI_Abort(MPI_COMM_WORLD, -99);
  } 
  
  if (me == 0 && !generateGraph && !isUnitEdgeWeight) {
      std::cerr << "Must specify -g for graph generation first before setting edge weights." << std::endl;
      MPI_Abort(MPI_COMM_WORLD, -99);
  }
  
  if (me == 0 && generateGraph && ((randomEdgePercent < 0) || (randomEdgePercent >= 100))) {
      std::cerr << "Invalid random edge percentage for generated graph!" << std::endl;
      MPI_Abort(MPI_COMM_WORLD, -99);
  }
} // parseCommandLine
