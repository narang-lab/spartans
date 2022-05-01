#include "TetLookup.h"
#include <cfloat>
#include <algorithm>
#include <cassert>

static const int nBins = 10;

TetLookup::TetLookup(const std::vector<vector3>& verts, std::vector<Tetrahedron>& tets) : verts(verts), tets(tets)
{
	//Determine min and max edge lengths:
	eMin = DBL_MAX; double eMax = 0.;
	for(const Tetrahedron& tet: tets)
		for(int j1=1; j1<4; j1++)
			for(int j2=0; j2<j1; j2++)
			{	double e = (verts[tet.iv[j1]] - verts[tet.iv[j2]]).length();
				if(e < eMin) eMin = e;
				if(e > eMax) eMax = e;
			}
	eMin -= 1e-6*(eMax-eMin); //so that eMin is a strict lower bound
	eMax += 1e-6*(eMax-eMin); //so that eMax is a strict upper bound
	deInv = nBins / (eMax - eMin);
}

//---------- O(1) tetrahedron comparison lookup -------------

template<int n> struct TetCompare; //Compare tetrahdera with n edges yet to be compared in length

template<int n> struct TetCompare //compare the n'th shortest edge
{	std::vector< TetCompare<n-1> > arr;
	
	//Find equivalent tetrahedron (also given sorted array of edge lengths, e in [0,nBins))
	const Tetrahedron* find(Tetrahedron& tet, const std::array<double,6>& e, const vector3* verts, double eTol, double vTolSq)
	{	if(!arr.size()) return 0;
		int eiMin = std::max(int(floor(e[6-n]-eTol)), 0);
		int eiMax = std::min(int(floor(e[6-n]+eTol)), nBins-1);
		for(int ei=eiMin; ei<=eiMax; ei++)
		{	const Tetrahedron* ret = arr[ei].find(tet, e, verts, eTol, vTolSq);
			if(ret) return ret;
		}
		return 0; //not found
	}
	
	//Add this tetrahedron for future lookups
	void add(const Tetrahedron& tet, const std::array<double,6>& e)
	{	if(!arr.size()) arr.resize(nBins);
		int ei = floor(e[6-n]);
		arr[ei].add(tet, e);
	}
};

template<> struct TetCompare<0> //recursion end: actually compare the vertex positions
{	std::vector< std::pair<std::array<double,6>, const Tetrahedron*> > cache;
	
	//Find equivalent tetrahedron (also given sorted array of edge lengths, e in [0,nBins))
	const Tetrahedron* find(Tetrahedron& tet, const std::array<double,6>& e, const vector3* verts, double eTol, double vTolSq)
	{	for(const auto& entry: cache)
		{	//First check edge lengths:
			const std::array<double,6>& ePrev = entry.first;
			bool match = true;
			for(int n=0; n<6; n++)
				if(fabs(e[n] - ePrev[n]) > eTol)
				{	match = false;
					break;
				}
			if(!match) continue;
			//Next check orientation (in all permutations):
			const Tetrahedron& tetPrev = *(entry.second);
			vector4i iv; //permuted vertices of tet
			for(int j0=0; j0<4; j0++)
			{	iv[0] = tet.iv[j0];
				vector3 xRef[4];
				xRef[0] = verts[iv[0]] - verts[tetPrev.iv[0]];
				for(int j=1; j<4; j++)
					xRef[j] = xRef[0] + verts[tetPrev.iv[j]];
				//Match all permutations of remaining vertices against xRef:
				for(int j1=0; j1<4; j1++) if(j1!=j0)
				{	iv[1] = tet.iv[j1];
					if((verts[iv[1]] - xRef[1]).length_squared() > vTolSq) continue;
					for(int j2=0; j2<4; j2++) if(j2!=j0 && j2!=j1)
					{	iv[2] = tet.iv[j2];
						if((verts[iv[2]] - xRef[2]).length_squared() > vTolSq) continue;
						int j3 = 6 - (j0+j1+j2);
						iv[3] = tet.iv[j3];
						if((verts[iv[3]] - xRef[3]).length_squared() > vTolSq) continue;
						//All vertex displacements match:
						tet.iv = iv; //updae the vertex permutation
						return &tetPrev;
					}
				}
			}
		}
		return 0; //not found
	}
	
	//Add this tetrahedron for future lookups
	void add(const Tetrahedron& tet, const std::array<double,6>& e)
	{	cache.push_back(std::make_pair(e,&tet));
	}
};

int TetLookup::assignIndices(double tol)
{	int nUnique = 0;
	TetCompare<6> tc; //compare 6 edge lengths and then actual tet geometry
	const double eTol = tol*nBins;
	const double vTolSq = std::pow(tol*(eMin+nBins/deInv), 2);
	std::vector< std::pair<size_t,size_t> > freq(tets.size()); //count frequency of tet repetition
	std::vector<bool> isUnique(tets.size(), false);
	for(size_t iTet=0; iTet<tets.size(); iTet++)
	{	Tetrahedron& tet = tets[iTet];
		//Calculate sorted edge length array:
		std::array<double,6> e; int ie=0;
		for(int j1=1; j1<4; j1++)
			for(int j2=0; j2<j1; j2++)
				e[ie++] = deInv*((verts[tet.iv[j1]] - verts[tet.iv[j2]]).length() - eMin);
		std::sort(e.begin(), e.end());
		//Search for equivalent tet:
		const Tetrahedron* tetEquiv = tc.find(tet, e, verts.data(), eTol, vTolSq);
		if(tetEquiv) //equivalent one already exists
			tet.iTetUnique = tetEquiv->iTetUnique;
		else //no previous equivalent
		{	tet.iTetUnique = iTet;
			freq[iTet].second = iTet;
			isUnique[iTet] = true;
			tc.add(tet, e); //add to list
			nUnique++;
		}
		freq[tet.iTetUnique].first++;
	}
	//Resort tetrahedra so that the prototypes of the most frequent ones appear first:
	std::sort(freq.begin(), freq.end(), std::greater< std::pair<size_t,size_t> >()); //descending order
	std::vector<Tetrahedron> tetsNew; tetsNew.reserve(tets.size());
	//--- add the unique ones first:
	for(int i=0; i<nUnique; i++)
	{	Tetrahedron& tet = tets[freq[i].second];
		tet.iTetUnique = i;
		tetsNew.push_back(tet);
	}
	//--- add the remaining ones:
	for(size_t iTet=0; iTet<tets.size(); iTet++)
		if(!isUnique[iTet])
		{	Tetrahedron& tet = tets[iTet];
			tet.iTetUnique = tets[tet.iTetUnique].iTetUnique;
			tetsNew.push_back(tet);
		}
	std::swap(tets, tetsNew);
	return nUnique;
}
