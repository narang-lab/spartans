#include "Util.h"
#include <map>

const StopWatch* stopWatchManager(const StopWatch* addWatch=0, const std::string* watchName=0)
{	static std::multimap<std::string, const StopWatch*> watches; //static array of all watches
	if(watchName)
	{	if(addWatch) watches.insert(std::make_pair(*watchName,addWatch));
		else //find and return watch pointer:
		{	auto iter = watches.find(*watchName);
			if(iter==watches.end()) return 0;
			else return iter->second;
		}
	}
	else //print timings:
	{	printf("\n");
		for(const auto& wPair: watches) wPair.second->print();
	}
	return 0;
}

StopWatch::StopWatch(std::string name) : Ttot(0), TsqTot(0), nT(0), name(name)
{	stopWatchManager(this, &name);
}

void StopWatch::start()
{	tPrev = clock_us();
}
void StopWatch::stop()
{	double T = clock_us()-tPrev;
	Ttot+=T; TsqTot+=T*T; nT++;
}

void StopWatch::print() const
{	if(nT)
	{	double meanT = Ttot/nT;
		double sigmaT = sqrt(TsqTot/nT - meanT*meanT);
		printf("PROFILER: %30s %12.6lf +/- %12.6lf s, %4d calls, %13.6lf s total\n",
			name.c_str(), meanT*1e-6, sigmaT*1e-6, nT, Ttot*1e-6);
	}
}

void StopWatch::printAll()
{	stopWatchManager();
}

void StopWatch::getTimings(std::string name, int& nT, double& Ttot)
{	const StopWatch* watch = stopWatchManager(0, &name);
	if(watch)
	{	nT = watch->nT;
		Ttot = watch->Ttot;
	}
	else
	{	nT = 0;
		Ttot = 0.;
	}
}
