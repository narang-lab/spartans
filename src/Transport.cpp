#include "OverlapInv.h"
#include "Util.h"
#include <iostream>
#include <string>
#include <cassert>
#include <algorithm>
#include <mpi.h>
#include <cfloat>

int transport(MPI_Comm comm, std::string inFile, std::string outFile)
{
	//MPI_Init(&argc, &argv);
	int nProcs, iProc;
	MPI_Comm_size(comm, &nProcs);
	MPI_Comm_rank(comm, &iProc);
	if(iProc>0) assert(freopen("/dev/null", "w", stdout));

	//Inputs:
	//StopWatch watchInput("Input"); watchInput.start();

	//--- open file:
	hid_t plid = H5Pcreate(H5P_FILE_ACCESS);
	H5Pset_fapl_mpio(plid, comm, MPI_INFO_NULL);
	hid_t fid = H5Fopen(inFile.c_str(), H5F_ACC_RDONLY, plid);
	H5Pclose(plid);
	if(fid<0) { fprintf(stderr, "Could not open input HDF5 file\n"); MPI_Abort(comm, 1); }

	//--- mesh:
	MeshTransport mt;
	mt.readInputs(comm, fid);
	mt.initialize(); //Initialize mesh

	//--- event list:
	std::vector<vector3> kArr;
	std::vector< std::vector<double> > BiArr;
	std::vector< std::vector<vector3> > SiArr;
	size_t nEvents = 0;
	int nBins = 0;

	printf("Explicit event mode.\n");
	kArr = h5readMatrix<3,double>(comm,fid,"k"); //decay wavevectors
	BiArr = h5readVectorArray<double>(comm, fid,"Bi");//body inputs
	SiArr = h5readMatrixArray<3,double>(comm, fid,"Si");//surface inputs

	//Check inputs:
	nEvents = kArr.size();
	nBins   = kArr.size();
	assert(BiArr.size() == size_t(nBins));
	for(const std::vector<double>& Bi: BiArr)
		if(Bi.size()) assert(Bi.size() == mt.verts.size());
	assert(SiArr.size() == size_t(nBins));
	for(const std::vector<vector3>& Si: SiArr)
		if(Si.size()) assert(Si.size() == mt.extSindex.size());

	H5Fclose(fid);
	//watchInput.stop();
	printf("nBins: %d  nEvents: %lu  nVerts: %lu\n", nBins, nEvents, mt.verts.size());

	//Initialize event lists
	//StopWatch watchEventPrep("EventListPrep"); watchEventPrep.start();

	//Split events evenly over MPI:
	size_t iEventStart = (nEvents * iProc) / nProcs;
	size_t iEventStop = (nEvents * (iProc+1)) / nProcs;
	std::vector<int> nEventsProc(nProcs, 0); nEventsProc[iProc] = iEventStop - iEventStart;
	std::vector<int> binEvents(nBins, 0); //number of events by energy bin
	std::vector<int> procStart(nBins, nProcs); //lowest process to involve bin location
	std::vector<int> procStop(nBins, -1); //highest process to involve bin location

	for(size_t iEvent=iEventStart; iEvent<iEventStop; iEvent++)
	{
		int ie = iEvent;
		binEvents[ie]++;
		//Mark all the bins that are active on this process:
		procStart[ie] = iProc;
		procStop[ie] = iProc;
	}

	MPI_Allreduce(MPI_IN_PLACE, nEventsProc.data(), nProcs, MPI_INT, MPI_SUM, comm);
	MPI_Allreduce(MPI_IN_PLACE, binEvents.data(), nBins, MPI_INT, MPI_SUM, comm);
	MPI_Allreduce(MPI_IN_PLACE, procStart.data(), nBins, MPI_INT, MPI_MIN, comm);
	MPI_Allreduce(MPI_IN_PLACE, procStop.data(), nBins, MPI_INT, MPI_MAX, comm);
	std::vector<int> nBinsProc(nProcs, 0);
	for(int iBin=0; iBin<nBins; iBin++)
	{	if(procStart[iBin]==nProcs) //not found on any proc i.e. no events in surrounding intervals
		{	assert(procStop[iBin]==-1);
			//Handle empty bin in last process:
			//procStart[iBin] = procStop[iBin] = (iBin ? procStop[iBin-1] : 0);
		}
		else assert(procStop[iBin]>=0);
		//update bin counts by process:
		for(int jProc=procStart[iBin]; jProc<=procStop[iBin]; jProc++)
			nBinsProc[jProc]++;
	}
	printf("nEvents by process: ["); for(int n: nEventsProc) printf(" %d", n); printf(" ]\n");
	//watchEventPrep.stop();

	//Make BiArr available where necessary according to this division:
	for(int ie=0; ie<nBins; ie++)
	{	std::vector<double>& Bi = BiArr[ie];
		std::vector<vector3>& Si = SiArr[ie];
		if(Bi.size()) //data owner
		{	for(int jProc=procStart[ie]; jProc<=procStop[ie]; jProc++)
				if(jProc != iProc)
				{	MPI_Send(Bi.data(), Bi.size(), MPI_DOUBLE, jProc, ie, comm);
					MPI_Send(Si.data(), Si.size()*3, MPI_DOUBLE, jProc, ie, comm);
				}
			if(!(procStart[ie]<=iProc && iProc<=procStop[ie]))
			{	Bi.clear(); //no longer needed on this process
				Si.clear();
			}
		}
		else if(procStart[ie]<=iProc && iProc<=procStop[ie]) //not data owner but need it
		{	Bi.resize(mt.verts.size());
			Si.resize(mt.extSindex.size());
			MPI_Recv(Bi.data(), Bi.size(), MPI_DOUBLE, MPI_ANY_SOURCE, ie, comm, MPI_STATUS_IGNORE);
			MPI_Recv(Si.data(), Si.size()*3, MPI_DOUBLE, MPI_ANY_SOURCE, ie, comm, MPI_STATUS_IGNORE);
		}
		//else don't have data and don't need it (do nothing)
	}
	

	//Prepare outputs;
	std::vector<double> intBSiArr(nBins, 0.);
	std::vector< std::vector<double> > OBoArr(nBins);
	std::vector< std::vector<vector3> > SoArr(nBins);
	std::vector< std::vector<double> > inOut(nBins);
	for(int ie=0; ie<nBins; ie++)
		if(procStart[ie]<=iProc && iProc<=procStop[ie])
		{	OBoArr[ie].resize(mt.verts.size());
			SoArr[ie].resize(mt.extSindex.size());
			inOut[ie].resize(mt.extSindex.size());
		}

	//Loop over events:
	int eventInterval = std::max(1, int(round(nEventsProc[iProc]/50.))); //interval for reporting progress
	int nEventsDone = 0; int nFailed = 0;
	
	printf("\nProcessing events: "); fflush(stdout);
	for(size_t iEvent=iEventStart; iEvent<iEventStop; iEvent++)
	{	//Figure out energy bin:
		int ie = iEvent;
		//Prepare input:
		const vector3& k = kArr[iEvent];
		std::vector<double> Bi(mt.verts.size());
		std::vector<vector3> Si(mt.extSindex.size());

		for(size_t iv=0; iv<mt.verts.size(); iv++)
			Bi[iv] = BiArr[ie][iv];

		// track carriers attempting to leave structure
		int flippedSignSi = 0;
		for(size_t iFace=0; iFace<mt.extSindex.size(); iFace++)
		{	double cosTheta = -dot(k,mt.extSnormals[iFace])/k.length(); //angle to inward normal
			//double wTheta = (cosTheta>0. ? 4*cosTheta : 0.); //contributions only on faces on input side
			
			double wTheta = 0.0;
			if(cosTheta>0.0) {
				wTheta = 1.0;
				// This is 1.0 and not cosTheta since both Si and So are stored 
				// per unit actual area, and mt.transport deals with
				// projecting/deprojecting.
			} else
			{	wTheta = 0.0;
				if (SiArr[ie][iFace][0] != 0.0 || SiArr[ie][iFace][1] != 0.0 || SiArr[ie][iFace][2] != 0.0) {
					flippedSignSi++;
				}
			}

			Si[iFace] = wTheta * SiArr[ie][iFace];
			inOut[ie][iFace] = wTheta;
		}

		//printf("\nWarning: Zeroed out input carriers attempting to exit structure, %d times.\n", flippedSignSi);
		
		double intSi = mt.integral(Si);
		double intBSi = mt.integral(Bi) + intSi;
		intBSiArr[ie] += intBSi;

		//Perform transport:
		std::vector<vector3> So; std::vector<double> OBo;
		if(!mt.transport(k, Si, Bi, So, OBo)) nFailed++;
		//Store output (with linear-spline histogramming):
		for(size_t iv=0; iv<mt.verts.size(); iv++)
		{	OBoArr[ie][iv] += OBo[iv];
		}
		for(size_t is=0; is<mt.extSindex.size(); is++)
		{	SoArr[ie][is] += So[is];
		}
		//Print progress:
		if((nEventsDone+1) % eventInterval == 0)
		{	printf("%d%% ", int(round((nEventsDone+1)*100./nEventsProc[iProc])));
			fflush(stdout);
		}
		nEventsDone++;
	}
	printf("done.\n"); fflush(stdout);

	//Communicate results over MPI:
	for(int ie=0; ie<nBins; ie++)
	{	//For each bin, the lowest involved process is responsible for final procesing and output:
		if(iProc==procStart[ie])
		{	//Recieve contributions from other processes:
			std::vector<double> OBoBuf(mt.verts.size());
			std::vector<vector3> SoBuf(mt.extSindex.size());
			std::vector<double> inOutBuf(mt.extSindex.size());
			for(int jProc=iProc+1; jProc<=procStop[ie]; jProc++)
			{	MPI_Recv(OBoBuf.data(), mt.verts.size(), MPI_DOUBLE, jProc, 0, comm, MPI_STATUS_IGNORE);
				MPI_Recv(SoBuf.data(), 3*mt.extSindex.size(), MPI_DOUBLE, jProc, 1, comm, MPI_STATUS_IGNORE);
				MPI_Recv(inOutBuf.data(), mt.extSindex.size(), MPI_DOUBLE, jProc, 1, comm, MPI_STATUS_IGNORE);
				for(size_t iv=0; iv<mt.verts.size(); iv++) OBoArr[ie][iv] += OBoBuf[iv];
				for(size_t is=0; is<mt.extSindex.size(); is++) 
				{	SoArr[ie][is] += SoBuf[is];
					inOut[ie][is] += inOutBuf[is];
				}
			}
		}
		if(procStart[ie]<iProc && iProc<=procStop[ie]) //send to responsible process
		{	MPI_Send(OBoArr[ie].data(), mt.verts.size(), MPI_DOUBLE, procStart[ie], 0, comm);
			MPI_Send(SoArr[ie].data(), 3*mt.extSindex.size(), MPI_DOUBLE, procStart[ie], 1, comm);
			MPI_Send(inOut[ie].data(), mt.extSindex.size(), MPI_DOUBLE, procStart[ie], 1, comm);
			OBoArr[ie].clear(); //keep data only on one proc
			SoArr[ie].clear(); //keep data only on one proc
			inOut[ie].clear(); //keep data only on one proc
		}
	}
	MPI_Allreduce(MPI_IN_PLACE, &nFailed, 1, MPI_INT, MPI_SUM, comm);
	if(nFailed) printf("WARNING: %d of %lu events skipped due to degeneracy issues.\n", nFailed, nEvents);

	//Final processing and output:
	OverlapInv oInv(mt);
	std::vector<double> intBoArr(nBins, 0.), intSoArr(nBins, 0.);
	printf("Applying Oinv: "); fflush(stdout);
	FILE* nullLog = fopen("/dev/null", "w");
	std::vector< std::vector<double> > BoArr(nBins);
	for(int ie=0; ie<nBins; ie++) if(iProc==procStart[ie])
	{	const std::vector<vector3>& So = SoArr[ie];
		std::vector<double>& Bo = BoArr[ie];
		Bo = oInv.apply(OBoArr[ie], 20, 1e-6, nullLog);
		//Save integrals for sum rule check:
		intBoArr[ie] = mt.integral(Bo);
		intSoArr[ie] = mt.integral(So);
	}
	fclose(nullLog);
	printf("done.\n");

	//Output:
	//StopWatch watchOutput("Output"); watchOutput.start();
	//--- open file:
	plid = H5Pcreate(H5P_FILE_ACCESS);
	H5Pset_fapl_mpio(plid, comm, MPI_INFO_NULL);
	fid = H5Fcreate(outFile.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, plid);
	if(fid<0) { fprintf(stderr, "Could not open/create output HDF5 file\n"); MPI_Abort(comm, 1); }
	H5Pclose(plid);
	//--- write data:
	h5writeVectorArray(comm, fid, "Bo", BoArr); //Bo data
	h5writeMatrixArray(comm, fid, "So", SoArr); //So data
	h5writeVectorArray(comm, fid, "inOut", inOut); //InOut data
	h5writeMatrix(comm, fid, "SoIndex", mt.extSindex); //vertex indices per face in So output
	//--- write diagnostic integrals:
	MPI_Allreduce(MPI_IN_PLACE, intBSiArr.data(), nBins, MPI_DOUBLE, MPI_SUM, comm);
	h5writeVector(comm, fid, "intBSi", intBSiArr); //Bi+Si; spatial integral
	//--- close file:
	H5Fclose(fid);
	//watchOutput.stop();

	//Sum rule check:
	printf("\nSum rule check:\n"); fflush(stdout);
	MPI_Allreduce(MPI_IN_PLACE, intSoArr.data(), nBins, MPI_DOUBLE, MPI_SUM, comm);
	MPI_Allreduce(MPI_IN_PLACE, intBoArr.data(), nBins, MPI_DOUBLE, MPI_SUM, comm);
	double largestSumRuleErr = 0.;
	int iLargestSumRuleErr = 0;
	for(int ie=0; ie<nBins; ie++)
	{	if(!intBSiArr[ie]) continue;
		double currentSumRuleErr = (intSoArr[ie]+intBoArr[ie])/intBSiArr[ie] - 1.;
	        //printf("IntSo: %le, IntBo: %le, IntBSi: %le, error: %le, at ie: %d\n", intSoArr[ie], intBoArr[ie], intBSiArr[ie], currentSumRuleErr, ie);
		if(currentSumRuleErr>largestSumRuleErr)
		{	largestSumRuleErr = currentSumRuleErr;
			iLargestSumRuleErr = ie;
		}
	}
	printf("Largest sum rule error: %le, at ie: %d\n", largestSumRuleErr,iLargestSumRuleErr);
	fflush(stdout);

	//Mtets report:
	//int nCalls; double Ttot; StopWatch::getTimings("Transport", nCalls, Ttot);
	//double Mtets = nCalls * mt.tets.size() / Ttot;
	//MPI_Allreduce(MPI_IN_PLACE, &nCalls, 1, MPI_INT, MPI_SUM, comm);
	//MPI_Allreduce(MPI_IN_PLACE, &Ttot, 1, MPI_DOUBLE, MPI_SUM, comm);
	//double MtetsAvg = nCalls * mt.tets.size() / Ttot;
	//double MtetsTot = MtetsAvg * nProcs;
	//printf("\nThroughput(Mtets): %.2lf local, %.2lf avg/process, %.2lf total.\n", Mtets, MtetsAvg, MtetsTot);

	//MPI_Finalize();
	//StopWatch::printAll();
	if(iProc>0) fclose(stdout);
	return 0;
}
