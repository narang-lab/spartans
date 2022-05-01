#include "MeshTransport.h"
#include "Util.h"
#include "TetLookup.h"
#include <cfloat>
#include <cassert>

void MeshTransport::readInputs(MPI_Comm comm, hid_t fid)
{	verts = h5readMatrix<3,double>(comm, fid, "verts");
	std::vector<vector4i> tets_iv = h5readMatrix<4,int>(comm, fid, "tets");
	tets.resize(tets_iv.size());
	for(size_t iTet=0; iTet<tets.size(); iTet++)
		tets[iTet].iv = tets_iv[iTet];
}

//Check for degeneracy in projected areas:
inline bool isProjectionDegenerate(const vector4i& iv, const std::vector<vector3>& verts, const double Atol)
{	//Get projected displacements of three vertices from the fourth:
	vector2 origin(verts[iv[3]]);
	vector2 dx[3];
	for(int j=0; j<3; j++)
		dx[j] = vector2(verts[iv[j]]) - origin;
	//Calculate areas of all four triangle:
	vector4 A;
	A[0] = dx[1][0]*dx[2][1] - dx[2][0]*dx[1][1];
	A[1] = dx[2][0]*dx[0][1] - dx[0][0]*dx[2][1];
	A[2] = dx[0][0]*dx[1][1] - dx[1][0]*dx[0][1];
	A[3] = A[0] + A[1] + A[2]; //working with signed areas so far
	double Asum = 0;
	for(int j=0; j<4; j++)
	{	A[j] = fabs(A[j]); //make unsigned
		Asum += A[j];
	}
	//Check areas relative to total projected:
	double Athresh = Asum * Atol;
	for(int j=0; j<4; j++)
		if(A[j]<Athresh)
			return true;
	return false;
}


//Get the face index containing specificed vertices (create and add this face to array, if not already in it)
inline int getFace(const vector3i& iv, std::vector<Face>& faces, std::vector< std::vector<int> >& vFaces)
{	// Find the vertex in iv which appears in the fewest known faces.
	// We pick this one to start with because it makes checking
	// for matching faces faster (fewer to check).
	size_t minFaces = vFaces[iv[0]].size(); int jMin=0;
	for(int j=1; j<3; j++)
	{	size_t nFaces = vFaces[iv[j]].size();
		if(nFaces < minFaces)
		{	minFaces = nFaces;
			jMin = j;
		}
	}
	//Check if the other two vertices are in any of the faces:
	for(int f: vFaces[iv[jMin]])
	{	const vector3i& fv = faces[f].iv;
		// Already know that iv[jMin] is present in fv,
		// because we picked faces that contained the vertex iv[jMin].
		// Check if the other two vertices in iv are in fv.
		bool matched = true;
		for(int j=0; j<3; j++) if(j!=jMin)
		{	bool found = false;
			for(int k=0; k<3; k++)
				if(iv[j] == fv[k])
				{	found = true;
					break;
				}
			if(!found)
			{	matched = false;
				break;
			}
		}
		if(matched) return f; //face already exists
	}
	//Create face:
	int fNew = faces.size(); //the index of the new face
	faces.push_back(Face(iv));
	for(int j=0; j<3; j++)
		vFaces[iv[j]].push_back(fNew); //update vertex->face map
	return fNew;
}


void MeshTransport::initialize()
{	//static StopWatch watchTetLookup("Init::tetLookup"), watchFaceArray("Init::faceArray");
	
	//Tetrahedron uniqueness:
	//watchTetLookup.start();
	int nUniqueTets = TetLookup(verts, tets).assignIndices();
	VtetUnique.resize(nUniqueTets);
	for(int iTet=0; iTet<nUniqueTets; iTet++)
	{	const Tetrahedron& tet = tets[iTet];
		vector3 x3 = verts[tet.iv[3]];
		matrix33 T;
		for(int j=0; j<3; j++)
			T.set_row(j, verts[tet.iv[j]] - x3);
		VtetUnique[iTet] = (1./6)*fabs(det(T));
	}
	printf("nTets: %lu  nUniqueTets: %d  expectedSpeedup: %.lfx\n", tets.size(), nUniqueTets, tets.size()*1./nUniqueTets);
	//watchTetLookup.stop();
	
	//Initialize global face array (with vertex maps):
	//watchFaceArray.start();
	faces.clear();
	faces.reserve(tets.size()*3); //conservative estimate
	std::vector< std::vector<int> > vFaces(verts.size(), std::vector<int>());
	std::vector<bool> isFaceExt(tets.size()*4, false); //absolute upper bound on size; bit array so memory not an issue
		
	// For the specified face in tet (defined in the scope of the next loop)
	// check to see if we've already encountered that face before in the loop
	// by directly examining the indices of the relevant vertices (in an order-agnostic
	// manner). If we haven't, add the face to the global face list.
	// Either way, store the global face index of that face in tet.iFace[f].
	#define SETFACE(f,j0,j1,j2) \
	{	int iFace = getFace(vector3i(tet.iv[j0],tet.iv[j1],tet.iv[j2]), faces, vFaces); \
		isFaceExt[iFace] = !isFaceExt[iFace]; \
		tet.iFace[f] = iFace; \
	}
	for(Tetrahedron& tet: tets) // Loop over tets
	{ // Process each face. The first index is the index of the face in the tetrahedron,
	  // the next three indices specify which vertices of the tetrahedron form that face.
		SETFACE(0, 1,2,3)
		SETFACE(1, 2,3,0)
		SETFACE(2, 3,0,1)
		SETFACE(3, 0,1,2)
	}
	#undef SETFACE
	vFaces.clear(); //no longer needed

	//--- external faces:
	size_t nExtFaces = 2*faces.size() - 4*tets.size();
	extSindex.clear();
	extSindex.reserve(nExtFaces);
	for(size_t iFace=0; iFace<faces.size(); iFace++)
		if(isFaceExt[iFace])
		{	Face& face = faces[iFace];
			face.iExt = extSindex.size();
			extSindex.push_back(face.iv);
		}
	isFaceExt.clear(); //no longer needed
	assert(extSindex.size() == nExtFaces);

	//--- external unit normals:
	extSnormals.resize(nExtFaces);
	for(const Tetrahedron& tet: tets)
		for(int f=0; f<4; f++)
		{	const Face& face = faces[tet.iFace[f]];
			if(face.iExt >= 0)
			{	//Get vector from interior vertex to first face vertex:
				int ivInt = dot(tet.iv,vector4i(1,1,1,1)) - dot(face.iv,vector3i(1,1,1)); //interior vertex index
				vector3 dOut = verts[face.iv[0]] - verts[ivInt]; //known to be an outward direction
				//Get normal direction:
				vector3 n = cross(verts[face.iv[0]]-verts[face.iv[1]], verts[face.iv[0]]-verts[face.iv[2]]);
				//Make sure n points outwards:
				if(dot(n,dOut) < 0.)
					n = -n; // Makes us agnostic to the ordering of the vertices in the face.
				extSnormals[face.iExt] = n / n.length(); //save normalized version
			}
		}

	//--- external face areas:
	extSareas.resize(nExtFaces);
	for(size_t f=0; f<extSindex.size(); f++)
	{	vector3 origin(verts[extSindex[f][2]]);
		vector3 dx0 = verts[extSindex[f][0]] - origin;
		vector3 dx1 = verts[extSindex[f][1]] - origin;
		extSareas[f] = 0.5*cross(dx0,dx1).length();
	}
	printf("nFaces: %lu  nExtFaces: %lu\n", faces.size(), nExtFaces);
	//watchFaceArray.stop();
}

bool MeshTransport::transport(vector3 k, const std::vector<vector3>& Si, const std::vector<double>& Bi, std::vector<vector3>& So, std::vector<double>& OBo) const
{	//static StopWatch watchTotal("Transport"), watchRotate("Transport::rotate"), watchSubTet("Transport::subTet"), watchMain("Transport::main");
	//watchTotal.start();
	
	//Rotate everything to align k along z, and if necessary,
	//joggle k a bit to prevent degenerate convex hull issues in subTet division
	//watchRotate.start();
	bool nonDegenerate = false;
	const int nJoggle = 100;
	const double Atol = 1e-8; //tolerance on relative areas to determine degeneracy
	const double dAngle = 1e-6; //typical rotation angle for degeneracy breaking attempts
	vector2 dkNet; //total tangential perturbation to k
	std::vector<vector3> verts = this->verts; //don't modify global vertex array (so that it can be reused)
	for(int iJoggle=0; iJoggle<nJoggle; iJoggle++)
	{
		//Calculate rotation that alignes decay vector to z:
		//--- extract norm and unit vector:
		double kNorm = k.length();
		vector3 kHat = (1./kNorm) * k;
		//--- construct an orthonormal set including kHat:
		int iAxis = -1; double kHatMinComp = DBL_MAX;
		for(int i=0; i<3; i++)
			if(fabs(kHat[i]) < kHatMinComp)
			{	kHatMinComp = fabs(kHat[i]);
				iAxis = i;
			}
		vector3 kPerp(0,0,0); kPerp[iAxis] = 1.; //axis with least overlap with kHat
		kPerp -= kHat * dot(kHat, kPerp); //orhogonalize w.r.t kHat
		kPerp *= (1./kPerp.length());
		//--- compute a rotation matrix that aligns k to the z direction:
		matrix33 rotInv;
		rotInv.set_col(0, cross(kPerp,kHat));
		rotInv.set_col(1, kPerp);
		rotInv.set_col(2, kHat);
		matrix33 rot = inv(rotInv);
		
		//Transform k and vertex array:
		k = rot * k;
		for(vector3& vert: verts)
			vert = rot * vert;
		
		//Check for degeneracy in the xy projection convex hull:
		nonDegenerate = true;
		for(const Tetrahedron& tet: tets)
			if(isProjectionDegenerate(tet.iv, verts, Atol))
			{	nonDegenerate = false;
				break;
			}
		if(nonDegenerate)
		{	if(iJoggle)
				printf("Joggled k by %.1le radians with %d moves to remove projected convex hull degeneracies\n", dkNet.length()/kNorm, iJoggle);
			break;
		}
		
		//Add small xy components to slightly rotate k:
		const double dkMag = kNorm * dAngle;
		for(int j=0; j<2; j++)
		{	k[j] = dkMag * Random::normal();
			dkNet[j] += k[j];
		}
	}
	//watchRotate.stop();
	
	//Split tetrahedra for current propagation direction and initialize tet queue:
	//watchSubTet.start();
	std::vector<int> check; check.reserve(3*tets.size()); //index of tetrahedra to check (which could be ready for processing)
	std::vector<bool> queued(tets.size(), false);
	std::vector<bool> ready(faces.size(), false);
	std::vector<int> iTetIn(faces.size(), -1);
	std::vector<int> iTetOut(faces.size(), -1);
	std::vector<SubTetArray> subTetsArr(VtetUnique.size());
	std::vector<vector4b> inFaceArr(VtetUnique.size()); //whether each face is incoming or outgoing on each tet
	std::vector<vector3> S(faces.size()); //intermediate values on the faces (used in main propagation loop)
	for(size_t iTet=0; iTet<VtetUnique.size(); iTet++)
		inFaceArr[iTet] = tets[iTet].setSubTets(k[2], verts.data(), subTetsArr[iTet]);
        for(size_t iTet=0; iTet<tets.size() && nonDegenerate; iTet++)
	{	const Tetrahedron& tet = tets[iTet];
		const vector4b& inFace = inFaceArr[tet.iTetUnique];
		//Update iTet's in the face array:
		for(int j=0; j<4; j++)
		{	int iFace = tet.iFace[j];
			int& iTetTarget = inFace[j] ? iTetOut[iFace] : iTetIn[iFace]; //note that this is the tetOut for an inFace and vice versa
			if(iTetTarget>=0) //iTet should be unique
			{	nonDegenerate = false;
				break;
			}
			iTetTarget = iTet;
			//Mark face ready if global input:
			const Face& face = faces[iFace];
			if(face.iExt>=0 && inFace[j])
			{	ready[iFace] = true;
				if(!queued[iTet]) //add tet to queue (if not already in it)
				{	check.push_back(iTet);
					queued[iTet] = true;
				}
				//Add surface input:
				double A = face.area(verts.data()); //actual area
				double Aproj = projectedArea(verts[face.iv[0]], verts[face.iv[1]], verts[face.iv[2]]);
				S[iFace] = Si[face.iExt] * A / fabs(Aproj); //convert from flux per actual area to flux per projected area
			}
		}
	}
	//watchSubTet.stop();
	
	//Main propagation loop:
	OBo.assign(verts.size(), 0.);
	So.assign(extSindex.size(), vector3());
	if(!nonDegenerate)
	{	printf("Failed to fix degeneracy after %d joggle attempts.\n", nJoggle);
		//watchTotal.stop();
		return false;
	}
	//watchMain.start();
	size_t nDone = 0;
	for(size_t iCheck=0; iCheck<check.size(); iCheck++)
	{	int iTet = check[iCheck];
		const Tetrahedron& tet = tets[iTet]; //current tetrahedon
		const vector4b& inFace = inFaceArr[tet.iTetUnique];
		//Check if all inputs are ready:
		bool allReady = true;
		for(int j=0; j<4; j++)
			if(inFace[j] && !ready[tet.iFace[j]])
			{	allReady = false;
				break;
			}
		if(allReady) //Propagate:
		{	//Collect surface inputs as 4-vectors in tet vertex order:
			vector4 Stet[4];
			for(int j=0; j<4; j++)
				if(inFace[j])
				{	int iFace = tet.iFace[j];
					const Face& face = faces[iFace];
					for(int fi=0; fi<3; fi++) //face vertex number
						for(int ti=0; ti<4; ti++) //tet vertex number
							if(face.iv[fi] == tet.iv[ti])
							{	Stet[j][ti] = S[iFace][fi];
								break;
							}
				}
			//Collect body inputs in tet vertex order:
			vector4 BiTet, BoTet;
			for(int j=0; j<4; j++)
				BiTet[j] = Bi[tet.iv[j]];
			//Loop over subtets:
			double Aout[4] = {0.,0.,0.,0.}; //projected areas of output faces
			const SubTetArray& subTets = subTetsArr[tet.iTetUnique];
			for(const SubTet& st: subTets)
			{	vector4 ziBot = zi3to4(st.ziBot);
				vector4 ziTop = zi3to4(st.ziTop);
				//Fetch the surface input:
				vector3 Si;
				const vector4& SiTet = Stet[st.fBot];
				Si[0] = SiTet[st.ia];
				Si[1] = SiTet[st.ib];
				Si[2] = dot(SiTet, ziBot);
				//Fetch the body input:
				vector4 Bi;
				Bi[0] = BiTet[st.ia];
				Bi[1] = BiTet[st.ib];
				Bi[2] = dot(BiTet, ziBot);
				Bi[3] = dot(BiTet, ziTop);
				//Propagate:
				vector3 So = st.A * (subTets.SiSo * Si + subTets.BiSo * Bi); // Here So is dimensionless, Si has units of 1/area. Later we apply invAfac to So to restore units of 1/area.
				vector4 Bo = st.A * (subTets.SiBo * Si + subTets.BiBo * Bi); // Here Bo is dimensionless, Bi has units of 1/volume. Later we apply the inverse overlap integral to Bo to restore units of 1/volume.
				//Surface output to tet vertex order:
				vector4& SoTet = Stet[st.fTop];
				SoTet[st.ia] += So[0];
				SoTet[st.ib] += So[1];
				SoTet += ziTop * So[2];
				//Body output to tet vertex order:
				BoTet[st.ia] += Bo[0];
				BoTet[st.ib] += Bo[1];
				BoTet += ziBot * Bo[2];
				BoTet += ziTop * Bo[3];
				//Collect areas:
				Aout[st.fTop] += st.A;
			}
			//Push surface outputs to global face array:
			for(int j=0; j<4; j++)
				if(!inFace[j])
				{	int iFace = tet.iFace[j];
					const Face& face = faces[iFace];
					double invAfac = 3./Aout[j]; //diagonal inverse overlap operator, restores units of 1/area.
					for(int fi=0; fi<3; fi++) //face vertex number
						for(int ti=0; ti<4; ti++) //tet vertex number
							if(face.iv[fi] == tet.iv[ti])
							{	S[iFace][fi] = invAfac * Stet[j][ti];
								break;
							}
					//mark output faces ready
					assert(!ready[iFace]);
					ready[iFace] = true;
					if(iTetOut[iFace] < 0) //global output
					{	assert(face.iExt >= 0);
						So[face.iExt] = S[iFace] * Aout[j]/face.area(verts.data()); //convert from flux per projected area to flux per actual area
					}
					else //Add next tet to check queue:
					{	if(!queued[iTetOut[iFace]])
						{	check.push_back(iTetOut[iFace]);
							queued[iTetOut[iFace]] = true;
						}
					}
				}
			//Push body output to global body array:
			for(int j=0; j<4; j++)
				OBo[tet.iv[j]] += BoTet[j];
			nDone++;
		}
		else //Not yet ready, add to end of list to try again later
		{	check.push_back(iTet);
		}
	}
	assert(nDone == tets.size()); //all tets must be processed
	//watchMain.stop();

	//watchTotal.stop();
	return true;
}

double MeshTransport::integral(const std::vector<double>& B) const
{	double result = 0.;
	for(const Tetrahedron& tet: tets)
	{	double Vfrac = 0.25*VtetUnique[tet.iTetUnique];
		double Bsum = 0.;
		for(int j=0; j<4; j++)
			Bsum += B[tet.iv[j]];
		result += Vfrac * Bsum;
	}
	return result;
}

double MeshTransport::integral(const std::vector<vector3>& S) const
{	double result = 0.;
	for(size_t f=0; f<extSindex.size(); f++)
	{	double Ssum = 0.;
		for(int j=0; j<3; j++)
			Ssum += S[f][j];
		result += extSareas[f] * (1./3) * Ssum;
	}
	return result;
}

