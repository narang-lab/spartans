#ifndef TETRAHEDRONTRANSPORT_TETLOOKUP_H
#define TETRAHEDRONTRANSPORT_TETLOOKUP_H

#include "Tetrahedron.h"

class TetLookup
{
public:
	TetLookup(const std::vector<vector3>& verts, std::vector<Tetrahedron>& tets);
	int assignIndices(double tol=1e-6); //assign tetrahedra iTet and iTetUnique, and return number of unique tetrahedra
private:
	const std::vector<vector3>& verts; //vertex array
	std::vector<Tetrahedron>& tets; //tetrahedron array
	double eMin, deInv; //min edge length and 1/de in the edge lookup table
};

#endif //TETRAHEDRONTRANSPORT_TETLOOKUP_H