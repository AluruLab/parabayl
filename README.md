# parabayl
Parallel Exact Bayesian Learning for small Bayesian networks (around. 30 variables).

## Algorithm
Repository contains souce code for the algorithm presented in the paper: 

Nikolova, O., Zola, J., & Aluru, S. (2013). Parallel globally optimal structure learning of Bayesian networks. *Journal of Parallel and Distributed Computing*, 73(8), 1039-1048.

Link to paper : https://doi.org/10.1016/j.jpdc.2013.04.001

## Requirements

This software is implemented using C++ and MPI-2 standards.

**Jam Make Redux** is required for building the sofware. It is available with most of the modern operating systems. 
If your OS does not provide "jam", you can compile from Jam source code available here: http://freetype.sourceforge.net/jam/ or from here: http://www.perforce.com/jam/jam.html.

**C++ compiler** – Any C++98 standard conforming compiler.

**MPI library** – Any MPI-2 compliant library. Tested libraries include MPICH2, OpenMPI, MVAPICH2.

**Boost Graph Library** – Version 1.49 or newer is required.

## Installation

The default architecture is mpich. The building process for the selected
architecture is initialized as:

    $ ./x-build <ARCH>

where ARCH is the selected architecture. For example:

    $ ./x-build bluegene

will build ParaBayL for IBM BlueGene. If the build process is successful a bin directory
with **parabayl_bluegene** executable in it will be generated. To finalize the installation
process, copy the binary to your preferred location. (The default building command is:

    $ ./x-build

and will generate the executable parabayl. )

When re-building, use the following command to properly clean previous projects prior to running
the x-build command:

    $ ./x-build distclean
    
    
## Running

To run ParaBayL use the command:

    $ mpiexec –np <NUMP> <PATH_TO_BIN>/parabayl <INPUTF> <OUTPUTF> <D> [<SYNC>]

where:

_NUMP_: number of processors.

_PATH_TO_BIN_: the path to the binary directory.

_INPUTF_: input file, *.exp* format.
The input is in *.exp* format and examples of this format are given in the data folder. 

_OUTPUTF_: output file, *.sif* format
The output is in standard .*sif* format
(http://wiki.cytoscape.org/Cytoscape_User_Manual/Network_Formats)
where each edge is represented by an entry and the first and second columns correspond
to parents and children, respectively, e.g.:

    A pp B 
implies A --> B in the output BN.

    * pp A 
implies that A does not have any parents.


_D_: in-degree bound (1 <= D <= n-1).

_SYNC_: optional parameter for number of iterations after which
communication is synchronized; best scaling results are achieved
with a value of 1 (equivalent to blocking communication) for
large data


## License

Our code is licensed under the Apache License 2.0 (see LICENSE). The licensing does not apply to the src/mpix and src/jaz folders. Please refer to the individual files in the folder for their licensing terms.
