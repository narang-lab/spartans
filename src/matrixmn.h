#ifndef TETRAHEDRONTRANSPORT_MATRIXMN_H
#define TETRAHEDRONTRANSPORT_MATRIXMN_H

#include "vectorn.h"

//Shorthands for commonly used sizes:
template<int m, int n, typename scalar=double> class matrixmn;
typedef matrixmn<2,2> matrix22;
typedef matrixmn<2,3> matrix23;
typedef matrixmn<2,4> matrix24;
typedef matrixmn<3,2> matrix32;
typedef matrixmn<3,3> matrix33;
typedef matrixmn<3,4> matrix34;
typedef matrixmn<4,2> matrix42;
typedef matrixmn<4,3> matrix43;
typedef matrixmn<4,4> matrix44;
typedef matrixmn<3,3,complex> matrix33z;

#define LOOPn(code) { for(int k=0; k<n; k++) { code } }
#define LOOPm(code) { for(int k=0; k<m; k++) { code } }
#define LOOPmn(code) { for(int i=0; i<m; i++) for(int j=0; j<n; j++) { code } }

template<int m, int n, typename scalar> class matrixmn
{
	scalar data[m][n];

public:
	//accessors:
	__hostanddev__ scalar& operator()(int i, int j) { return data[i][j]; }
	__hostanddev__ const scalar& operator()(int i, int j) const { return data[i][j]; }
	__hostanddev__ vectorn<n,scalar> row(int i) const { vectorn<n,scalar> ret; LOOPn( ret[k]=data[i][k]; ) return ret; }
	__hostanddev__ vectorn<m,scalar> column(int j) const { vectorn<m,scalar> ret; LOOPm( ret[k]=data[k][j]; ) return ret; }

	__hostanddev__ void set_row(int i, const vectorn<n,scalar>& v) { LOOPn( data[i][k] = v[k]; ) }
	__hostanddev__ void set_col(int j, const vectorn<m,scalar>& v) { LOOPm( data[k][j] = v[k]; ) }

	//constructors:
	__hostanddev__ matrixmn(){ LOOPmn( data[i][j]=0; ) }
	template<typename scalar2> explicit __hostanddev__ matrixmn(const matrixmn<m,n,scalar2>& other) { LOOPmn( data[i][j] = scalar(other(i,j)); ) }
	
	//arithmetic operators
	__hostanddev__ matrixmn operator-() const { matrixmn ret; LOOPmn( ret(i,j) = -data[i][j]; ) return ret; }
	__hostanddev__ matrixmn operator+(const matrixmn &other) const { return (matrixmn(*this) += other); }
	__hostanddev__ matrixmn& operator+=(const matrixmn &other) { LOOPmn( data[i][j] += other(i,j); ) return *this; }
	__hostanddev__ matrixmn operator-(const matrixmn &other) const { return (matrixmn(*this) -= other); }
	__hostanddev__ matrixmn& operator-=(const matrixmn &other) { LOOPmn( data[i][j] -= other(i,j); ) return *this; }
	__hostanddev__ matrixmn& operator*=(scalar s) {	LOOPmn( data[i][j] *= s; ) return *this; }
	__hostanddev__ matrixmn operator*(scalar s) const { return (matrixmn(*this) *= s); }
	__hostanddev__ matrixmn& operator/=(scalar s) { return (*this) *= (1./s); }
	__hostanddev__ matrixmn operator/(scalar s) const { return (*this) * (1./s); }

	//! transpose
	__hostanddev__ matrixmn<n,m,scalar> operator~() const
	{	matrixmn<n,m,scalar> ret;
		LOOPmn( ret(j,i) = data[i][j]; )
		return ret;
	}

	void print(FILE* fp, const char *format, bool brackets=true) const
	{	for(int i=0; i<m; i++)
		{	if(brackets) fprintf(fp, "[ ");
			for(int j=0; j<n; j++) fprintf(fp, format, data[i][j]);
			if(brackets) fprintf(fp, " ]\n"); else fprintf(fp, "\n");
		}
	}

	//Comparison operators
	__hostanddev__ bool operator==(const matrixmn& other) const
	{	LOOPmn( if(data[i][j] != other(i,j)) return false; )
		return true;
	}
	__hostanddev__ bool operator!=(const matrixmn& other) const
	{	return ! (*this == other);
	}
};

//Multiplies:
template<int m, int n, typename scalar> __hostanddev__ matrixmn<m,n,scalar> operator*(scalar s, const matrixmn<m,n,scalar> &mat) { return mat*s; }

template<int m, int n, typename scalar> __hostanddev__ matrixmn<m,n,scalar> outer(const vectorn<m,scalar> &a, const vectorn<n,scalar> &b)
{	matrixmn<m,n,scalar> ret;
	LOOPmn( ret(i,j) = a[i] * b[j]; )
	return ret;
}

template<int m, int n, typename scalar> __hostanddev__ vectorn<m,scalar> operator*(const matrixmn<m,n,scalar>& mat, const vectorn<n,scalar> &v)
{	vectorn<m,scalar> ret;
	LOOPmn( ret[i] += mat(i,j) * v[j]; )
	return ret;
}

template<int m, int n, typename scalar> __hostanddev__ vectorn<n,scalar> operator*(const vectorn<m,scalar> &v, const matrixmn<m,n,scalar>& mat)
{	vectorn<n,scalar> ret;
	LOOPmn( ret[j] += v[i] * mat(i,j); )
	return ret;
}

template<int m, int k, int n, typename scalar> __hostanddev__ matrixmn<m,n,scalar> operator*(const matrixmn<m,k,scalar> &a, const matrixmn<k,n,scalar>& b)
{	matrixmn<m,n,scalar> ret;
	LOOPmn( for(int l=0; l<k; l++) ret(i,j) += a(i,l) * b(l,j); )
	return ret;
}

template<int m, int n, typename scalar> __hostanddev__ matrixmn<m,n,scalar>& operator*=(matrixmn<m,n,scalar> &mat, const matrixmn<n,n,scalar>& other)
{	return (mat = mat * other);
}


template<int n, typename scalar> __hostanddev__ matrixmn<n,n,scalar> Diag(vectorn<n,scalar> v)
{	matrixmn<n,n,scalar> ret;
	LOOPn( ret(k,k) = v[k]; )
	return ret;
}

template<int n, typename scalar> __hostanddev__ scalar trace(const matrixmn<n,n,scalar> &mat) { scalar ret=0; LOOPn( ret += mat(k,k); ) return ret; }
template<int m, int n> __hostanddev__ double nrm2(const matrixmn<m,n>& mat) { double ret=0; LOOPmn( ret += mat(i,j)*mat(i,j); ) return sqrt(ret); } //!< 2-norm of matrix
template<typename scalar> __hostanddev__ scalar det(const matrixmn<3,3,scalar> &mat) { return box(mat.row(0),mat.row(1),mat.row(2)); }
template<typename scalar> __hostanddev__ matrixmn<3,3,scalar> adjugate(const matrixmn<3,3,scalar> &mat)
{	matrixmn<3,3,scalar> adj;
	adj.set_col(0, cross(mat.row(1),mat.row(2)));
	adj.set_col(1, cross(mat.row(2),mat.row(0)));
	adj.set_col(2, cross(mat.row(0),mat.row(1)));
	return adj;
}
__hostanddev__ matrixmn<3,3> inv(const matrixmn<3,3> &mat)
{	return (1./det(mat)) * adjugate(mat);
}

#undef LOOPm
#undef LOOPn

#endif // TETRAHEDRONTRANSPORT_MATRIXMN_H
