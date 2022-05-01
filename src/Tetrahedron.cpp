#include "Tetrahedron.h"
#include <cfloat>
#include <cassert>
#include <iostream>

//-------- Misc. Helper functions -------

//Comparison checking for smallest in y, and then x
inline static bool yxLess(const vector3& a, const vector3& b)
{	return (a[1]<b[1]) || ( (a[1]==b[1]) && (a[0]<b[0]) );
}

//Get baryentric coordinates of a vertex:
inline static vector3 ziVertex(int i)
{	vector3 zi; //all zeroes, which corresponds to vertex i=3
	if(i<3) zi[i] = 1.;
	return zi;
}

//Convert face barycentric coordinates to tetrahedron ones:
inline static vector3 ziFacePoint(int i0, int i1, int i2, double t1, double t2)
{	vector3 zi;
	if(i0<3) zi[i0] = 1.-(t1+t2);
	if(i1<3) zi[i1] = t1;
	if(i2<3) zi[i2] = t2;
	return zi;
}

//Convert edge barycentric coordinates to tetrahedron ones:
inline static vector3 ziEdgePoint(int i0, int i1, double t)
{	vector3 zi;
	if(i0<3) zi[i0] = 1.-t;
	if(i1<3) zi[i1] = t;
	return zi;
}

//--------------- class Tetrahedron ----------

vector4b Tetrahedron::setSubTets(double k, const vector3* vArr, SubTetArray& subTets) const
{
	vector3 x[4];
	matrix33 T; //transformation from barycentric to Cartesian coordinates (without the overall shift)

	//Load vertices:
	for(int j=0; j<4; j++)
		x[j] = vArr[iv[j]];
	
	//Set transformation:
	for(int j=0; j<3; j++)
		T.set_col(j, x[j]-x[3]);

	//Calculate xy-plane convex hull (to distinguish cases):
	std::vector<int> iHull; iHull.reserve(4);
	std::vector<bool> onHull(4, false);
	//--- find bottom-most point:
	int iyMin = 0;
	for(int i=1; i<4; i++)
		if(yxLess(x[i], x[iyMin]))
			iyMin = i;
	iHull.push_back(iyMin); onHull[iyMin] = true;
	//--- add remaining points by Jarvis march algorithm:
	vector3 dPrev(1.,0.,0.); //direction of previous segment
	while(true)
	{	int iCur = iHull.back();
		int iNext = -1;
		double cosNext = -DBL_MAX, dotNext = DBL_MAX;
		for(int i=0; i<4; i++)
			if(i != iCur)
			{	vector3 dx = x[i] - x[iCur];
				double dot = dPrev[0]*dx[0] + dPrev[1]*dx[1];
				double cross = dPrev[0]*dx[1] - dPrev[1]*dx[0];
				double cos = dot / hypot(cross, dot);
				//Select the point that has lowest angle to previous segment (and shortest distance if same angle)
				if(cos>cosNext || (cos==cosNext && dot<dotNext))
				{	iNext = i;
					cosNext = cos;
					dotNext = dot;
				}
			}
		if(iNext == iyMin) break; //convex hull has been closed
		iHull.push_back(iNext); onHull[iNext] = true;
		dPrev = x[iNext] - x[iCur];
	}
	assert(iHull.size()==3 || iHull.size()==4);
	
	//Split into triangular projections based on convex hull:
	subTets.resize(iHull.size());
	if(iHull.size()==3)
	{	//Projection contains one vertex "tip" inside a triangle
		//--- get tip:
		int iTip = -1;
		for(int i=0; i<4; i++)
			if(!onHull[i])
			{	iTip = i;
				break;
			}
		vector3 ziTip = ziVertex(iTip);
		//--- find point to split opposite face at:
		vector3 dx1 = x[iHull[1]] - x[iHull[0]];
		vector3 dx2 = x[iHull[2]] - x[iHull[0]];
		vector3 dxTip = x[iTip] - x[iHull[0]];
		//--- --- solve dx1 t1 + dx2 t2 = dxTip projected to xy-plane:
		double den = 1./(dx1[0]*dx2[1] - dx1[1]*dx2[0]);
		double t1 = den * (dxTip[0]*dx2[1] - dxTip[1]*dx2[0]);
		double t2 = den * (dx1[0]*dxTip[1] - dx1[1]*dxTip[0]);
		vector3 ziBase = ziFacePoint(iHull[0], iHull[1], iHull[2], t1, t2);
		//--- figure out input and output directions:
		double zTip = x[iTip][2];
		double zBase = x[iHull[0]][2] + t1*dx1[2] + t2*dx2[2];
		bool tipBot = (zTip < zBase);
		subTets.h = fabs(zTip - zBase);
		//--- add a SubTet for each edge of convex hull:
		for(int e=0; e<3; e++)
		{	SubTet& st = subTets[e];
			st.ia = iHull[e];
			st.ib = iHull[e==2 ? 0 : e+1];
			int fTip = iHull[e ? e-1 : 2]; //index of face containing tip, which is opposite the vertex left out in base
			const int& fBase = iTip; //since base is on face opposite tip
			if(tipBot)
			{	st.fBot = fTip;
				st.fTop = fBase;
				st.ziBot = ziTip;
				st.ziTop = ziBase;
			}
			else
			{	st.fBot = fBase;
				st.fTop = fTip;
				st.ziBot = ziBase;
				st.ziTop = ziTip;
			}
			st.A = projectedArea(x[iTip], x[st.ia], x[st.ib]);
		}
	}
	else //iHull.size()==4
	{	//Projection is a quadrilateral
		//--- calculate intersection of xy-projected diagonals
		//'even' diagonal: x[iHull[0]] (1-tE) + tE x[iHull[2]]
		//'odd' diagonal:  x[iHull[1]] (1-tO) + tO x[iHull[3]]
		//therefore solve xy projection of tE dxE - tO dxO = dx10 with the definitions:
		vector3 dxE = x[iHull[2]] - x[iHull[0]];
		vector3 dxO = x[iHull[3]] - x[iHull[1]];
		vector3 dx10 = x[iHull[1]] - x[iHull[0]];
		double den = 1./(dxE[0]*dxO[1] - dxE[1]*dxO[0]);
		double tE =  den * (dx10[0]*dxO[1] - dx10[1]*dxO[0]);
		double tO = -den * (dxE[0]*dx10[1] - dxE[1]*dx10[0]);
		//--- get coordinates of intersection points:
		vector3 ziE = ziEdgePoint(iHull[0], iHull[2], tE);
		vector3 ziO = ziEdgePoint(iHull[1], iHull[3], tO);
		vector3 xE = x[iHull[0]] + tE*dxE;
		double zE = xE[2];
		double zO = x[iHull[1]][2] + tO*dxO[2]; //xy components of xO not needed as they are same as xE
		bool oddBot = (zO < zE);
		subTets.h = fabs(zO - zE);
		//--- add a SubTet for each edge of convex hull:
		for(int e=0; e<4; e++)
		{	SubTet& st = subTets[e];
			st.ia = iHull[e];
			st.ib = iHull[e==3 ? 0 : e+1];
			int fE = 6 - (st.ia + st.ib + iHull[(e==1|| e==2) ? 0 : 2]); //face index containing even diagonal and current edge
			int fO = 6 - (st.ia + st.ib + iHull[(e==0|| e==1) ? 3 : 1]); //face index containing odd diagonal and current edge
			if(oddBot)
			{	st.fBot = fO;
				st.fTop = fE;
				st.ziBot = ziO;
				st.ziTop = ziE;
			}
			else
			{	st.fBot = fE;
				st.fTop = fO;
				st.ziBot = ziE;
				st.ziTop = ziO;
			}
			st.A = projectedArea(xE, x[st.ia], x[st.ib]);
		}
	}
	
	//checkSubTets(subTets, x, T); //uncomment to test the result of above code
	subTets.transportInit(k); //initialize transport coefficients
	
	//Update in/out of faces:
	vector4b inFace;
	for(const SubTet& st: subTets)
	{	inFace[st.fBot] = true;
		inFace[st.fTop] = false;
	}
	return inFace;
}

void Tetrahedron::checkSubTets(const SubTetArray& subTets, const vector3 (&x)[4], const matrix33& T) const
{	double vol = (1./6) * fabs(det(T));
	std::cout << "vol: " << vol << "  nHull: " << subTets.size() << "  checking subTets";
	const double tol = 1e-8;
	double volSum = 0.;
	for(size_t i=0; i<subTets.size(); i++)
	{	std::cout << ' ' << i; std::cout.flush();
		const SubTet& st = subTets.at(i);
		//Check that the faces include the common edge:
		assert(st.ia != st.fTop); assert(st.ia != st.fBot);
		assert(st.ib != st.fTop); assert(st.ib != st.fBot);
		//Check that the top and bot vertices lie on specified faces:
		if(st.fTop<3) assert(fabs(st.ziTop[st.fTop]) < tol); else assert(fabs(1.-dot(st.ziTop,vector3(1,1,1))) < tol);
		if(st.fBot<3) assert(fabs(st.ziBot[st.fBot]) < tol); else assert(fabs(1.-dot(st.ziBot,vector3(1,1,1))) < tol);
		//Check that the top vertex is directy above bot vertex:
		vector3 xTop = x[3] + T * st.ziTop;
		vector3 xBot = x[3] + T * st.ziBot;
		assert(fabs(xTop[0] - xBot[0]) < tol);
		assert(fabs(xTop[1] - xBot[1]) < tol);
		assert(xTop[2] > xBot[2]);
		//Verify area and height:
		assert(fabs(subTets.h - (xTop[2]-xBot[2])) < tol);
		assert(fabs(st.A - projectedArea(x[st.ia], x[st.ib], xTop)) < tol);
		//Increment sum for volume check:
		volSum += (1./3) * st.A * subTets.h;
	}
	double volSumRelErr = fabs(volSum/vol - 1.);
	std::cout << "  volSumRelErr: " << volSumRelErr << '\n'; std::cout.flush();
	assert(volSumRelErr < tol);
}
