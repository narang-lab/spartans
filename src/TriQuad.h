#ifndef TETRAHEDRONTRANSPORT_TRIQUAD_H
#define TETRAHEDRONTRANSPORT_TRIQUAD_H

#include "vectorn.h"
#include <array>

//Entry in triangle quadrature:
struct TriQuadEntry
{	double w; //weight
	vector3 zi; //barycentric coordinates
	int nPerm; //number of permutations in symmetrization
	
	TriQuadEntry(double w) : w(w), zi(1./3,1./3,1./3), nPerm(1) {} //centroid node
	TriQuadEntry(double w, double a) : w(w), zi(a,a,1.-2.*a), nPerm(3) {} //median node
	TriQuadEntry(double w, double a, double b) : w(w), zi(a,b,1.-(a+b)), nPerm(6) {} //general node
};


extern std::array<vector3i,6> permutations; //Permutations of (0,1,2)

extern std::array<TriQuadEntry,1> triQuad3; //3-point quadrature correct to order 2
extern std::array<TriQuadEntry,2> triQuad6; //6-point quadrature correct to order 4
extern std::array<TriQuadEntry,3> triQuad7; //7-point quadrature correct to order 5
extern std::array<TriQuadEntry,3> triQuad12; //12-point quadrature correct to order 6

#endif // TETRAHEDRONTRANSPORT_TRIQUAD_H
