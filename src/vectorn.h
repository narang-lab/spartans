#ifndef TETRAHEDRONTRANSPORT_VECTORN_H
#define TETRAHEDRONTRANSPORT_VECTORN_H

#include <vector>
#include <cstdio>
#include <cmath>
#include <complex>
typedef std::complex<double> complex; //NOTE: this disallows "using namespace std" which I never do anyway

#define LOOPn(code) { for(int k=0; k<n; k++) { code } }
#define __hostanddev__ inline //potentialy add __device__ later on for GPU support

//Shorthands for commonly used sizes:
template<int n, typename scalar=double> class vectorn;
typedef vectorn<2> vector2;
typedef vectorn<3> vector3;
typedef vectorn<4> vector4;
typedef vectorn<3,complex> vector3z;
typedef vectorn<2,int> vector2i;
typedef vectorn<3,int> vector3i;
typedef vectorn<4,int> vector4i;
typedef vectorn<2,bool> vector2b;
typedef vectorn<3,bool> vector3b;
typedef vectorn<4,bool> vector4b;

template<int n, typename scalar> class vectorn
{
	scalar data[n];
public:
	//Accessors
	__hostanddev__ scalar& operator[](int k) { return data[k]; }
	__hostanddev__ const scalar& operator[](int k) const { return data[k]; }
	__hostanddev__ vectorn operator[](const vectorn<n,int>& idx) const { vectorn ret; LOOPn( ret[k] = data[idx[k]]; ) return ret; }
	
	//Constructor
	__hostanddev__ vectorn() { LOOPn(data[k]=0;) }
	vectorn(const std::vector<scalar>& a) { LOOPn(data[k]=a[k];) }
	template<int n2, typename scalar2> __hostanddev__ explicit vectorn(const vectorn<n2,scalar2>& a) { LOOPn(data[k]=(k<n2 ? a[k] : 0);) } //conversion between scalar types and lengths (fill unknown with zero)
	
	//Variable argument constructor using template recursion:
	inline void construct_helper(int k) {}
	template<typename... Args> inline void construct_helper(int k, scalar dk, Args... args) { data[k]=dk; construct_helper(k+1, args...); }
	template<typename... Args> inline vectorn(Args... args) { construct_helper(0, args...); }
	
	//Arithmetic:
	__hostanddev__ vectorn operator+(const vectorn &a) const { return (vectorn(*this) += a); }
	__hostanddev__ vectorn operator+=(const vectorn &a) { LOOPn( data[k]+=a[k]; ) return *this; }
	__hostanddev__ vectorn operator+(const scalar a) const { return (vectorn(*this) += a); }
	__hostanddev__ vectorn operator+=(const scalar a) { LOOPn( data[k]+=a; ) return *this; }

	__hostanddev__ vectorn operator-() const { vectorn ret(*this); LOOPn( ret[k]=-ret[k]; ) return ret; }
	__hostanddev__ vectorn operator-(const vectorn &a) const { return (vectorn(*this) -= a); }
	__hostanddev__ vectorn operator-=(const vectorn &a) { LOOPn( data[k]-=a[k]; ) return *this; }

	__hostanddev__ vectorn operator/(scalar s) const { return (*this)*(1.0/s); }
	__hostanddev__ vectorn& operator/=(scalar s) { return (*this)*=(1.0/s); }

	__hostanddev__ scalar length_squared() const { scalar ret=0; LOOPn( ret += data[k]*data[k]; ) return ret; }
	__hostanddev__ scalar length() const { return sqrt(length_squared()); }

	void print(FILE* fp, const char *format) const { std::fprintf(fp, "[ "); LOOPn( fprintf(fp, format, data[k]); ) std::fprintf(fp, " ]\n"); }
	__hostanddev__ bool operator==(const vectorn& w) const { LOOPn( if(data[k] != w[k]) return false; ) return true; }
	__hostanddev__ bool operator<(const vectorn& w) const { LOOPn( if(data[k]!=w[k]) return data[k]<w[k]; ) return false; }
};

template<int n, typename scalar> __hostanddev__ vectorn<n,scalar> operator+(scalar s, const vectorn<n,scalar>& a) { return a+s; }
template<int n, typename scalar> __hostanddev__ vectorn<n,scalar>& operator*=(vectorn<n,scalar>& a, scalar s) { LOOPn(a[k]*=s;) return a; }
template<int n, typename scalar> __hostanddev__  vectorn<n,scalar> operator*(scalar s, const vectorn<n,scalar> &a) { vectorn<n,scalar> v; LOOPn(v[k]=a[k]*s;) return v; }
template<int n, typename scalar> __hostanddev__  vectorn<n,scalar> operator*(const vectorn<n,scalar> &a, scalar s) { vectorn<n,scalar> v; LOOPn(v[k]=a[k]*s;) return v; }
template<int n, typename scalar> __hostanddev__ scalar dot(const vectorn<n,scalar>& a, const vectorn<n,scalar>& b) { scalar ret=0; LOOPn( ret += a[k]*b[k]; ) return ret; }
template<int n> __hostanddev__ vectorn<n,complex> conj(const vectorn<n,complex>& a) { vectorn<n,complex> v; LOOPn(v[k]=conj(a[k]);) return v; }

//! cross product
template<typename scalar> __hostanddev__ vectorn<3,scalar> cross(const vectorn<3,scalar> &a, const vectorn<3,scalar> &b)
{	return vectorn<3,scalar>(
		a[1]*b[2] - a[2]*b[1],
		a[2]*b[0] - a[0]*b[2],
		a[0]*b[1] - a[1]*b[0] );
}

//! box product / triple product
template<typename scalar> __hostanddev__ scalar box(const vectorn<3,scalar>& a, const vectorn<3,scalar>& b, const vectorn<3,scalar>& c)
{	return dot(a,cross(b,c));
}

#undef LOOPn

#endif // TETRAHEDRONTRANSPORT_VECTORN_H
