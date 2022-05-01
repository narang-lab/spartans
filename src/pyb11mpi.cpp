#include "Transport.h"

#include <mpi4py/mpi4py.h>
#include <pybind11/pybind11.h>

#include <string>

namespace py = pybind11;

MPI_Comm *get_mpi_comm(py::object py_comm) {
  auto comm_ptr = PyMPIComm_Get(py_comm.ptr());

  if (!comm_ptr)
    throw py::error_already_set();

  return comm_ptr;
}

PYBIND11_MODULE(_core, m) {
  if (import_mpi4py() < 0) {
    throw py::error_already_set();
  }
  m.def(
      "transport",
      [](py::object py_comm, std::string inFile, std::string outFile) {
        auto comm = get_mpi_comm(py_comm);
        transport(*comm, inFile, outFile);
      },
      R"pbdoc(
           Perform transport.
    )pbdoc");
}
