#ifndef TETRAHEDRONTRANSPORT_OVERLAPINV_H
#define TETRAHEDRONTRANSPORT_OVERLAPINV_H

#include "MeshTransport.h"

class OverlapInv
{
public:
	OverlapInv(const MeshTransport&);
	
	std::vector<double> apply(const std::vector<double>& in, int nIterations=20, double tol=1e-6, FILE* fpLog=stdout) const;
	std::vector<double> precondition(const std::vector<double>& in) const; //diagonal approximation to apply
private:
	int nVerts;
	std::vector<double> oInvDiag; //diagonal entries of preconditioner
	const std::vector<Tetrahedron>& tets; //tetrahedra in mesh
	const std::vector<double>& VtetUnique; //tetrahedron volumes
	std::vector<double> hessian(const std::vector<double>& in) const; //overlap operator (apply() is inverse of this)
};

#endif // TETRAHEDRONTRANSPORT_OVERLAPINV_H
