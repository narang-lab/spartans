#include "OverlapInv.h"
#include "Util.h"

OverlapInv::OverlapInv(const MeshTransport& mt) : nVerts(mt.verts.size()), oInvDiag(nVerts, 0.), tets(mt.tets), VtetUnique(mt.VtetUnique)
{
	//Initialize diagonal approx to O:
	for(const Tetrahedron& tet: tets)
	{	const double Vfrac = 0.25 * VtetUnique[tet.iTetUnique];
		for(int j=0; j<4; j++)
			oInvDiag[tet.iv[j]] += Vfrac;
	}
	//Invert it to get preconditioner:
	for(double& x: oInvDiag)
		x = 1./x;
}

inline std::vector<double>& operator*=(std::vector<double>& X, double s) { for(double& x: X) x *= s; return X; }
inline void axpy(double alpha, const std::vector<double>& X, std::vector<double>& Y) { for(size_t i=0; i<X.size(); i++) Y[i] += alpha * X[i]; }
inline double dot(const std::vector<double>& X, const std::vector<double>& Y) { double ret=0.; for(size_t i=0; i<X.size(); i++) ret += X[i]*Y[i]; return ret; }

std::vector<double> OverlapInv::apply(const std::vector<double>& in, int nIterations, double tol, FILE* fpLog) const
{	//static StopWatch watch("OverlapInv"); watch.start();
	
	#define sync
	//Initialize:
	std::vector<double> x = precondition(in); //since it is the approximate inverse
	double rdotz0 = sync(dot(x, in)); //tolerance is relative to this
	std::vector<double> r = in; axpy(-1.0, hessian(x), r); //residual r = rhs - A.x;
	std::vector<double> z = precondition(r), d = r; //the preconditioned residual and search direction
	double beta=0.0, rdotzPrev=0.0, rdotz = sync(dot(r, z));

	//Check initial residual
	double rzNorm = rdotz0 ? sqrt(fabs(rdotz)/rdotz0) : 0.;
	fprintf(fpLog, "OinvCG: Initial:  sqrt(|r.z|): %12.6le\n", rzNorm); fflush(fpLog);
	//if(rzNorm<tol) { fprintf(fpLog, "OinvCG: Converged sqrt(r.z)<%le\n", tol); fflush(fpLog); watch.stop(); return x; }
	if(rzNorm<tol) { fprintf(fpLog, "OinvCG: Converged sqrt(r.z)<%le\n", tol); fflush(fpLog); return x; }

	//Main loop:
	int iter;
	for(iter=0; iter<nIterations; iter++)
	{	//Update search direction:
		if(rdotzPrev)
		{	beta = rdotz/rdotzPrev;
			d *= beta; axpy(1.0, z, d); // d = z + beta*d
		}
		else d = z; //fresh search direction (along gradient)
		//Step:
		std::vector<double> w = hessian(d);
		double alpha = rdotz/sync(dot(w,d));
		axpy(alpha, d, x);
		axpy(-alpha, w, r);
		z = precondition(r);
		rdotzPrev = rdotz;
		rdotz = sync(dot(r, z));
		//Print info:
		double rzNorm = sqrt(fabs(rdotz)/rdotz0);
		fprintf(fpLog, "OinvCG: Iter: %3d  sqrt(|r.z|): %12.6le  alpha: %12.6le  beta: %13.6le\n",
			iter, rzNorm, alpha, beta); fflush(fpLog);
		//Check convergence:
		//if(rzNorm<tol) { fprintf(fpLog, "OinvCG: Converged sqrt(r.z)<%le\n", tol); fflush(fpLog); watch.stop(); return x; }
		if(rzNorm<tol) { fprintf(fpLog, "OinvCG: Converged sqrt(r.z)<%le\n", tol); fflush(fpLog); return x; }
	}
	#undef sync
	fprintf(fpLog, "OinvCG: Gradient did not converge within threshold in %d iterations\n", iter); fflush(fpLog);
	//watch.stop(); return x;
	return x;
}

std::vector<double> OverlapInv::precondition(const std::vector<double>& in) const
{	std::vector<double> out(nVerts);
	for(int i=0; i<nVerts; i++)
		out[i] = oInvDiag[i] * in[i];
	return out;
}


std::vector<double> OverlapInv::hessian(const std::vector<double>& in) const
{	std::vector<double> out(nVerts, 0.);
	for(const Tetrahedron& tet: tets)
	{	double scaleFac = 0.05 * VtetUnique[tet.iTetUnique];
		double contribSum = 0.;
		for(int j=0; j<4; j++)
		{	double contrib = scaleFac * in[tet.iv[j]];
			out[tet.iv[j]] += contrib;
			contribSum += contrib;
		}
		for(int j=0; j<4; j++)
		out[tet.iv[j]] += contribSum;
	}
	return out;
}
