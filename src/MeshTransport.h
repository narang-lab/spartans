#ifndef TETRAHEDRONTRANSPORT_MESHTRANSPORT_H
#define TETRAHEDRONTRANSPORT_MESHTRANSPORT_H

#include "Tetrahedron.h"
#include "H5io.h"
#include <istream>

class MeshTransport
{
public:
	//Input quantities: set before calling initialize
	std::vector<vector3> verts; //global vertex array
	std::vector<Tetrahedron> tets; //tetrahedra (with the vertex indices iv set at input)
	
	void readInputs(MPI_Comm comm, hid_t fid); //read inputs from HDF5 file (datasets called verts and tets)
	
	void initialize(); //set up faces and mesh adjacency
	std::vector<vector3i> extSindex; //vertex indices of external triangles (order for surface output)
	std::vector<vector3> extSnormals; //outward unit normals of external faces
	std::vector<double> extSareas; //areas of exterior faces
	std::vector<double> VtetUnique; //volumes of the unique tetrahedra
	
	//Propagate input (Bi) to outputs (Bo,So):
	bool transport(vector3 k, const std::vector<vector3>& Si, const std::vector<double>& Bi, std::vector<vector3>& So, std::vector<double>& OBo) const; 

	double integral(const std::vector<double>& B) const; //calculate integral of a scalar field on the vertices
	double integral(const std::vector<vector3>& S) const; //calculate integral of a flux field on the faces
private:
	std::vector<Face> faces; //global array of faces
};

#endif // TETRAHEDRONTRANSPORT_MESHTRANSPORT_H
