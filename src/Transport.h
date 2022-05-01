#pragma once
#include <mpi.h>
#include <string>

int transport(MPI_Comm comm, std::string inFile, std::string outFile);

