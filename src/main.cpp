/***
 *  $Id: main.cpp 566 2010-04-01 21:32:04Z olia $
 **
 *  File: main.cpp
 *  Created: May 22, 2009
 *
 *  Author: Jaroslaw Zola <jaroslaw.zola@gmail.com>
 *          Olga Nikolova <olga.nikolova@gmail.com>
 */

#include <jaz/math_add.hpp>
#include <mpix/MPI_env.hpp>

#include <cstdlib>
#include <iostream>
#include <limits>

#include "AppConfig.hpp"
#include "InputData.hpp"
#include "data_transform.hpp"
#include "parabayl.hpp"
#include "read_exp.hpp"

#include <fstream>

#ifdef BLUE_GENE
#  include <rts.h>
#endif // BLUE_GENE


bool check_cpu(const mpix::MPI_env& mpi_env) {
    const int MIN_CPU = 2;
    if (mpi_env.size() < MIN_CPU) return false;
    double kd = log2(mpi_env.size());
    unsigned int ki = static_cast<unsigned int>(kd);
    if ((kd - ki) > std::numeric_limits<double>::epsilon()) return false;
    return true;
} // check_cpu


MPI_Comm topo_create(MPI_Comm in) {
    int nnodes = 0;
    MPI_Comm_size(MPI_COMM_WORLD, &nnodes);

    int k = static_cast<int>(log2(nnodes));

    // each node has degree k
    int* index = new int[nnodes];
    for (int i = 0; i < nnodes; ++i) index[i] = (i + 1) * k;

    int* edges = new int[nnodes * k];

    unsigned int pos = 0;
    unsigned int mask = 0;

    for (int i = 0; i < nnodes; ++i) {
	mask = (1 << k);
	for (unsigned int j = 0; j < k; ++j) {
	    mask >>= 1;
	    unsigned int p = (i | mask);
	    if (i != p) edges[pos++] = p;
	}
	mask = (1 << k);
	for (unsigned int j = 0; j < k; ++j) {
	    mask >>= 1;
	    unsigned int p = (i & (~mask));
	    if (i != p) edges[pos++] = p;
	}
    } // for

    MPI_Comm comm;
    MPI_Graph_create(in, nnodes, index, edges, 1, &comm);

    delete[] edges;
    delete[] index;

    return comm;
} // topo_create



bool run(const mpix::MPI_env& mpi_env,
	 const std::string& infile, const std::string& outfile, int bound, long int step) {
    AppConfig app_conf;

    app_conf.input = infile;
    app_conf.output = outfile;
    app_conf.deg_bound = bound;
    app_conf.sync_step = step;

    app_conf.mem_report = false;
    // app_conf.mem_report = true;

    InputData data;
    bool res = read_exp(mpi_env, app_conf, data);
    if (res == false) return false;

    // check if we can handle this problem
    /*if ((sizeof(bit_t) * 8 - 1) <  data.g_nrow) {
	std::cout << master << error << "too many variables" << std::endl;
	return false;
    }*/

    if (data.g_nrow < log2(mpi_env.size())) {
	std::cout << master << error << "too many processors" << std::endl;
	return false;
    }

    // check if bound makes sense
    if (bound < 1) {
	std::cout << master << error << "d needs to be >= 1" << std::endl;
	return false;
    }

    std::cout << master << "* transforming data" << std::endl;
    
    standardize(data);
    
    // *** PRINT STANDARDIZED DATA ***
    /*
    std::stringstream ss;
    ss.clear();
    ss.str(""); 
    ss << app_conf.output.c_str() << ".stddata";
  
    std::ofstream f;
    f.open(ss.str().c_str());
    if(!f) std::cout << error << "standardized data file could not be created" << std::endl << std::flush; 
    
    for (unsigned int i = 0; i < data.g_nrow; ++i) {
      f << master << i << ":\t" << std::flush;
      double* tab = data.exp_table.row_begin(i);
      for (unsigned int j = 0; j < data.g_ncol; ++j) {
	f << master << static_cast<double>(tab[j]) << " " << std::flush;
      }
      f << master << std::endl << std::flush;
    }
    f.close();
    */
    std_discretize(data);

    // *** PRINT DISCRETIZED DATA ***
    /*
    std::stringstream ss;
    ss.clear();
    ss.str(""); 
    ss << app_conf.output.c_str() << ".stddata";
  
    std::ofstream f;
   
    ss.clear();
    ss.str("");
    ss << app_conf.output.c_str() << ".dscdata";

    f.clear();
    f.open(ss.str().c_str());
    if(!f) std::cout << error << "discretized data file could not be created" << std::endl << std::flush; 

    for (unsigned int i = 0; i < data.g_nrow; ++i) {
      f << master << i << ":\t" << std::flush;
      unsigned char* tab = data.D.begin() + i * data.g_ncol;
      for (unsigned int j = 0; j < data.g_ncol; ++j) {
	f << master << static_cast<int>(tab[j]) << " " << std::flush;
      }
      f << master << std::endl << std::flush;
    }
    f.close();
    */

    // free original data
    data.exp_table = jaz::numeric_table<double>();
    if (app_conf.mem_report == true) MEM_REPORT;

    // sync to have nice startup
    MPI_Barrier(mpi_env.comm());

    res = parabayl(mpi_env, app_conf, data);

    return res;
} // run


int main(int argc, char* argv[]) {
    std::ios_base::sync_with_stdio(false);
    std::cin.tie(0);

    MPI_Init(&argc, &argv);

    // here we need to map hypercube to 3D mesh
    // ...

    // create MPI environment
    // int rank;
    // MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    // if (rank == 0) topo_create(MPI_COMM_WORLD);
    // return MPI_Finalize();

    mpix::MPI_env mpi_env(MPI_COMM_WORLD);
    mpi_env.root(0);

    mpi_iom.mpi_env = &mpi_env;

    // check processors
    if (check_cpu(mpi_env) == false) {
	std::cout << master << error << "incorrect number of processors" << std::endl;
	return MPI_Finalize();
    }

    // check arguments
    if ((argc != 4) && (argc != 5)) {
	std::cout << master << "Usage: " << argv[0] << " infile outfile indegbound [syncstep]" << std::endl;
	return MPI_Finalize();
    }

    // say hello
    if (mpi_env.am_I_root() == true) {
	std::cout << FULL_NAME << "\n"
		  << SHORT_NAME << " " << MAJOR_VERSION << "." << MINOR_VERSION
		  << " coded 2009-2010 by Olga Nikolova and Jaroslaw Zola" << std::endl;

	std::cout << "ARGV:";
	for (int i = 1; i < argc; ++i) std::cout << " " << argv[i];
	std::cout << ", CPU: " << mpi_env.size();
	std::cout << master << ", bit_t: " << sizeof(bit_t) * 8 << " bits";
	std::cout << std::endl;
    }

    // let's run
    TIMER_START;

    int bound = 64;
    if (argc > 3) {
      bound = std::atol(argv[3]);
    }

    long int step = 0;
    if (argc == 5) {
      step = std::atol(argv[4]);
    }

    bool res = run(mpi_env, argv[1], argv[2], bound, step);

    if (res == false) {
	std::cout << master << error << "something went wrong!" << std::endl;
    } else {
	MPI_Barrier(mpi_env.comm());
	std::cout << master << "\nDone: " << timer(TIMER_GET) << std::endl;
    }

    return MPI_Finalize();
} // main
