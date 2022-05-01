#ifndef TETRAHEDRONTRANSPORT_SUBTET_H
#define TETRAHEDRONTRANSPORT_SUBTET_H

#include "matrixmn.h"
#include "TriQuad.h"

#define triQuad triQuad12 //select the quadrature

//Tetrahedral cells after splitting original tetrahedron by z projection
struct SubTet
{	int ia, ib; //indices of vertices of common edge (0-3; barycentric coord is implicit)
	int fTop, fBot; //tetrahedron face numbers that the top and bottom faces are part of
	vector3 ziTop, ziBot; //barycentric coordinates of third vertex on top and bottom faces
	double A; //projected triangle area
};

//A collection of SubTets from the same tetrahedron
struct SubTetArray : public std::vector<SubTet>
{	double h; //maximum height (= real distance between ziTop and ziBot) common to SubTets from same tetrahedron

	//Coordinate systems for the transport matrices:
	//   S terms: triangle barycentric coordinates in order a, b, bot/top
	//   B terms: SubTet barycentric coordinates in order a, b, bot, top
	matrix33 SiSo;
	matrix43 SiBo;
	matrix34 BiSo;
	matrix44 BiBo;
	void transportInit(double k); //calculate transport coefficients
	void transportCheck(double k); //check transport coefficients against brute force Monte Carlo calculation
};

//Barycentric coordinate conversion
inline vector4 zi3to4(vector3 ziIn)
{	vector4 ziOut(ziIn); //copy first three coordinates
	ziOut[3] = 1. - (ziOut[0] + ziOut[1] + ziOut[2]);
	return ziOut;
}

#endif //TETRAHEDRONTRANSPORT_SUBTET_H
