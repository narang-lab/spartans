#include "SubTet.h"
#include "Util.h"

void SubTetArray::transportInit(double k)
{	const double a = k*h; // Tet height times decay wavenumber.
	SiSo = matrix33();
	SiBo = matrix43();
	BiSo = matrix34();
	BiBo = matrix44();
	// triQuad is a quadrature over a tetrahedron
	for(const TriQuadEntry tqe: triQuad) // tqe is a point in the quadrature
	{	const double dA = tqe.w; //weight in the quadrature
		//Calculate necessary exponential function combinations:
		vector3 exp0_base, exp1_base, exp2_base, exp3_base, exp4_base;
		for(int j=0; j<3; j++)
		{	double y = a*tqe.zi[j]; // zi is the barycentric coordinate of the quadrature point
			exp0_base[j] = exp(-y);
			if(fabs(y) > 0.15)
			{	exp1_base[j] = (1. - exp0_base[j])/y;
				exp2_base[j] = (1. - exp1_base[j])/y;
				exp3_base[j] = (0.5 - exp2_base[j])/y;
				exp4_base[j] = (1./6 - exp3_base[j])/y;
			}
			else //Use a series expansion that carefully minimizes roundoff/truncation errors:
			{	exp4_base[j] = (1./24)*(1. - (1./5)*y*(1. - (1./6)*y*(1. - (1./7)*y*(1. - (1./8)*y*(1. - (1./9)*y*(1. - (1./10)*y))))));
				exp3_base[j] = 1./6 - y * exp4_base[j];
				exp2_base[j] = 0.5 - y * exp3_base[j];
				exp1_base[j] = 1. - y * exp2_base[j];
			}
		}
		for(int iPerm=0; iPerm<tqe.nPerm; iPerm++)
		{	const vector3i& p = permutations[iPerm];
			//Fetch relevant quantities after permutation
			vector3 zi = tqe.zi[p];
			double exp0 = exp0_base[p[2]];
			double exp1 = exp1_base[p[2]];
			double exp2 = exp2_base[p[2]];
			double exp3 = exp3_base[p[2]];
			double exp4 = exp4_base[p[2]];
			//SiSo:
			SiSo += (dA * exp0) * outer(zi,zi); // Dimensionless
			//SiBo:
			double hCur = h*zi[2];
			vector4 expIntBo( zi[0]*exp1, zi[1]*exp1, zi[2]*exp2, zi[2]*(exp1 - exp2) );
			SiBo += (dA * hCur * k) * outer(expIntBo, zi); // Dimensionless
			//BiSo:
			vector4 expIntBi( expIntBo[0], expIntBo[1], expIntBo[3], expIntBo[2] ); //SiBo and BiSo related by simultaneous t, z reversal
			BiSo += (dA * hCur) * outer(zi, expIntBi); // Units of length
			//BiBo:
			vector4 ziExt( zi[0], zi[1], zi[2], zi[2] );
			matrix44 expIntBiBo = outer(ziExt, ziExt);
			expIntBiBo(0,0) *= exp2;
			expIntBiBo(0,1) *= exp2;
			expIntBiBo(1,0) *= exp2;
			expIntBiBo(1,1) *= exp2;
			expIntBiBo(0,2) *= (exp2-exp3);
			expIntBiBo(1,2) *= (exp2-exp3);
			expIntBiBo(3,0) *= (exp2-exp3);
			expIntBiBo(3,1) *= (exp2-exp3);
			expIntBiBo(0,3) *= exp3;
			expIntBiBo(1,3) *= exp3;
			expIntBiBo(2,0) *= exp3;
			expIntBiBo(2,1) *= exp3;
			expIntBiBo(2,2) *= (exp3-exp4);
			expIntBiBo(3,3) *= (exp3-exp4);
			expIntBiBo(2,3) *= exp4;
			expIntBiBo(3,2) *= (exp2-2.*exp3+exp4);
			BiBo += (dA * hCur*hCur * k) * expIntBiBo; // Units of length
		}
	}
}


void SubTetArray::transportCheck(double k)
{	const double a = k*h;
	const int N = 100000; //number of MC points
	double dA = (1./(N*permutations.size()));
	matrix33 SiSoMC;
	matrix43 SiBoMC;
	matrix34 BiSoMC;
	matrix44 BiBoMC;
	for(int i=0; i<N; i++)
	{	//Pick a random point in the triangle:
		vector3 zi_base;
		zi_base[0] = Random::uniform();
		zi_base[1] = Random::uniform();
		zi_base[2] = 1.-(zi_base[0]+zi_base[1]);
		if(zi_base[2] < 0.) //if in the wrong half of unit square
		{	//Reflect along diagonal:
			zi_base[0] = 1. - zi_base[0];
			zi_base[1] = 1. - zi_base[1];
			zi_base[2] = -zi_base[2];
		}
		//Pick random z fractions:
		double zFrac = Random::uniform();
		double zFracIn = Random::uniform(0., 1.-zFrac); //z fraction at input for BiBo, with zFrac then corresponding to zFracOut - zFracIn (to reuse exponential below)
		double zFracOut = zFracIn + zFrac;
		//Pre-calculate exponentials:
		vector3 exp0_base, exp0z_base;
		for(int k=0; k<3; k++)
		{	exp0_base[k] = exp(-a*zi_base[k]);
			exp0z_base[k] = exp(-a*zi_base[k]*zFrac);
		}
		//Accumulate contribution from all permutations of this point:
		for(const vector3i& p: permutations)
		{	//Fetch quantities after permutation:
			vector3 zi = zi_base[p];
			double exp0 = exp0_base[p[2]];
			double exp0z = exp0z_base[p[2]];
			//SiSo:
			SiSoMC += outer(zi,zi) * (dA * exp0);
			//SiBo:
			double wSiBo = zi[2]*a; // = zi[2]*h * k. First part for length of segment in MC integral, second for rate of deposition of carriers
			vector4 ziBo(zi[0], zi[1], zi[2]*(1.-zFrac), zi[2]*zFrac);
			SiBoMC += outer(ziBo,zi) * (dA * exp0z * wSiBo);
			//BiSo:
			double wBiSo = zi[2]*h; //like wBiSo, but without the k
			vector4 ziBi(zi[0], zi[1], zi[2]*zFrac, zi[2]*(1.-zFrac)); //body point at 1-zFrac along altitude so that transport fraction is still zFrac
			BiSoMC += outer(zi,ziBi) * (dA * exp0z * wBiSo);
			//BiBo:
			double wBiBo = k * std::pow(zi[2]*h,2) * (1.-zFrac); //integration weight for the double integral over z, for the chosen MC scheme of drawing zFrac
			vector4 ziBin(zi[0], zi[1], zi[2]*(1.-zFracIn), zi[2]*zFracIn);
			vector4 ziBout(zi[0], zi[1], zi[2]*(1.-zFracOut), zi[2]*zFracOut);
			BiBoMC += outer(ziBout,ziBin) * (dA * exp0z * wBiBo);
		}
	}
	transportInit(k);
	printf("a: %9.3le  SiSoErr: %9.3le  SiBoErr: %9.3le  BiSoErr: %9.3le  BiBoErr: %9.3le\n", k*h,
		nrm2(SiSo-SiSoMC)/nrm2(SiSoMC), nrm2(SiBo-SiBoMC)/nrm2(SiBoMC),
		nrm2(BiSo-BiSoMC)/nrm2(BiSoMC), nrm2(BiBo-BiBoMC)/nrm2(BiBoMC) );
}
