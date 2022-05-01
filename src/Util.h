#ifndef TETRAHEDRONTRANSPORT_UTIL_H
#define TETRAHEDRONTRANSPORT_UTIL_H

#include <random>
#include <sys/time.h>
#include <thread>

namespace Random
{	inline std::mt19937_64& generator() { static std::mt19937_64 g; return g; }
	inline double uniform() { static std::uniform_real_distribution<double> uniformDist; return uniformDist(generator()); }
	inline double uniform(double start, double stop) { return start + (stop-start)*uniform(); }
	inline double normal() { static std::normal_distribution<double> normdist; return normdist(generator()); }
}

inline double clock_us()
{	timeval tv;
	gettimeofday(&tv,NULL);
	return ((tv.tv_sec & 0x1fffff) * 1e6) + tv.tv_usec;
}

class StopWatch
{
public:
	StopWatch(std::string name);
	void start();
	void stop();
	void print() const;
	static void printAll();
	static void getTimings(std::string name, int& nT, double& Ttot);
private:
	double tPrev, Ttot, TsqTot; int nT;
	std::string name;
};

//Templated wrapper for std::thread
template<typename Callable,typename ... Args> void threadLaunch(int nThreads, Callable* func, Args... args)
{	std::thread** tArr = new std::thread*[nThreads-1];
	for(int t=0; t<nThreads; t++)
	{	if(t<nThreads-1) tArr[t] = new std::thread(func, t, nThreads, args...);
		else (*func)(t, nThreads, args...);
	}
	for(int t=0; t<nThreads-1; t++)
	{	tArr[t]->join();
		delete tArr[t];
	}
	delete[] tArr;
}

#endif // TETRAHEDRONTRANSPORT_UTIL_H
