#ifndef TETRAHEDRONTRANSPORT_TETRAHEDRON_H
#define TETRAHEDRONTRANSPORT_TETRAHEDRON_H

#include "SubTet.h"

//Triangular face in input mesh, containing connectivity information and surface transport values
struct Face
{	vector3i iv; //vertex indices
	int iExt; //index into external face array (-1 if an internal face)
	
	Face(const vector3i& iv = vector3i(-1,-1,-1)) : iv(iv), iExt(-1) {}
	
	//Return unsigned area of face given vertex array:
	inline double area(const vector3* verts) const
	{	vector3 v0 = verts[iv[0]];
		return 0.5*cross(verts[iv[1]]-v0, verts[iv[2]]-v0).length();
	}
};


//Tetrahedron in input mesh
struct Tetrahedron
{	vector4i iv; //indices into a vertex array
	vector4i iFace; //indices of faces (opp each vertex in same order) into global face array
	int iTetUnique; //index of this tetrahedron in list of unique tetrahedra
	
	//call after setting iv and iFace, only on the unique tets
	//returns which faces are in on the tetrahedron
	vector4b setSubTets(double k, const vector3* vArr, SubTetArray& subTets) const;
private:
	void checkSubTets(const SubTetArray& subTets, const vector3 (&x)[4], const matrix33& T) const; //Verify the subTet split
};

//Calculated z-projected (signed) area:
inline double projectedArea(const vector3& x0, const vector3& x1, const vector3& x2)
{	vector3 dx1 = x1 - x0;
	vector3 dx2 = x2 - x0;
	return 0.5*(dx1[0]*dx2[1] - dx1[1]*dx2[0]);
}

#endif //TETRAHEDRONTRANSPORT_TETRAHEDRON_H
