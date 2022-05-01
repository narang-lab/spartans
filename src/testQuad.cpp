#include "SubTet.h"

//Average a function with signature "double f(vector3 x, args...)" over a triangle, where x are barycentric coordinates
//(Multiply by triangle area to get integral over triangle)
template<typename Func, typename... Args> double triAverage(const Func& f, Args... args)
{	double result = 0.;
	for(const TriQuadEntry& tqe: triQuad)
		for(int iPerm=0; iPerm<tqe.nPerm; iPerm++)
			result += tqe.w * f(tqe.zi[permutations[iPerm]], args...);
	return result;
}

//exponential of each barycentric coordinate
inline double exp_bary(const vector3& x, double a, int dir)
{	return exp(-a*x[dir]);
}

//analytical average of exp_bary on triangle (independent of dir by symmetry)
inline double expAvg(double a)
{	return fabs(a)<0.06
		? 1. - (1./3)*a*(1. - (1./4)*a*(1. - (1./5)*a*(1. - (1./6)*a*(1. - (1./7)*a*(1. - (1./8)*a)))))
		: (exp(-a) + a - 1.) / (0.5*a*a);
}

int main()
{
	for(double a = 0.015; a<10.; a*=1.05)
	{	double exact = expAvg(a);
		printf("%9.3le  %9.3le ", a, exact);
		for(int dir=0; dir<3; dir++)
		{	double approx = triAverage(exp_bary, a, dir);
			printf("  %+9.2le", approx/exact-1.);
		}
		printf("\n");
	}
	return 0;
}
