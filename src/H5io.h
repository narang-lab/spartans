#ifndef TETRAHEDRONTRANSPORT_H5IO_H
#define TETRAHEDRONTRANSPORT_H5IO_H

#include <hdf5.h>
#include "matrixmn.h"
#include <vector>
#include <cassert>
#include <cstring>

template<typename T> T h5readScalar(MPI_Comm comm, hid_t fid, const char* aname); //Read scalar from an attribute
template<typename T> std::vector<T> h5readVector(MPI_Comm comm, hid_t fid, const char* dname); //Read 1D array from a dataset
template<int n, typename T> std::vector< vectorn<n,T> > h5readMatrix(MPI_Comm comm, hid_t fid, const char* dname);  //Read 2D array from a dataset
template<int m, int n, typename T> std::vector< matrixmn<m,n,T> > h5readMatrix(MPI_Comm comm, hid_t fid, const char* dname);  //Read 3D array from a dataset
template<typename T> std::vector< std::vector<T> > h5readVectorArray(MPI_Comm comm, hid_t fid, const char* dname);  //Collective read array of arrays (net 2D) from a dataset (data distributed uniformly amongst processes)
template<int n, typename T> std::vector<std::vector<vectorn<n,T> > > h5readMatrixArray(MPI_Comm comm, hid_t fid, const char* dname); //Colectively read array of matrices (net 3D) from a dataset (data distributed uniformly amongst processes)

template<typename T> void h5writeVector(MPI_Comm comm, hid_t fid, const char* dname, const std::vector<T>& data); //Collectively write 1D array to a dataset whanl all the data is available on all the processes
template<int n, typename T> void h5writeMatrix(MPI_Comm comm, hid_t fid, const char* dname, const std::vector< vectorn<n,T> >& data); //Collectively write 2D array when all the data is available on all processes
template<typename T> void h5writeVectorArray(MPI_Comm comm, hid_t fid, const char* dname, const std::vector<std::vector<T> >& data); //Collectively write array of arrays (net 2D) where each process stores mutually exclusive outermost indices
template<int n, typename T> void h5writeMatrixArray(MPI_Comm comm, hid_t fid, const char* dname, const std::vector<std::vector<vectorn<n,T> > >& data); //Collectively write array of matrices (net 3D) where each process stores mutually exclusive outermost indices

//---------- Implementations

template<typename T> struct h5type;
template<> struct h5type<int> { static hid_t get() { return H5T_NATIVE_INT; } };
template<> struct h5type<double> { static hid_t get() { return H5T_NATIVE_DOUBLE; } };
template<> struct h5type<complex>
{	static hid_t get()
	{	static hid_t result = -1;
		if(result < 0)
		{	result = H5Tcreate (H5T_COMPOUND, sizeof(complex));
			H5Tinsert(result, "r", 0*sizeof(double), H5T_NATIVE_DOUBLE);
			H5Tinsert(result, "i", 1*sizeof(double), H5T_NATIVE_DOUBLE);
		}
		return result;
	}
};

template<typename T> T h5readScalar(MPI_Comm comm, hid_t fid, const char* aname)
{	T result;
	hid_t aid = H5Aopen_name(fid, aname);
	if(aid<0){ fprintf(stderr, "Could not read attribute '%s' from input file.\n", aname); MPI_Abort(comm, 1); }
	H5Aread(aid, h5type<T>::get(), &result);
	H5Aclose(aid);
	return result;
}

template<typename T> std::vector<T> h5readVector(MPI_Comm comm, hid_t fid, const char* dname)
{	hid_t did = H5Dopen1(fid, dname);
	if(did<0){ fprintf(stderr, "Could not read dataset '%s' from input file.\n", dname); MPI_Abort(comm, 1); }
	//Get data dimensions:
	hid_t dspace = H5Dget_space(did);
	int nDims = H5Sget_simple_extent_ndims(dspace);
	assert(nDims==1 || nDims==2);
	hsize_t dims[2];
	H5Sget_simple_extent_dims(dspace, dims, NULL);
	for(int iDim=1; iDim<nDims; iDim++)
		assert(dims[iDim] == 1);
	//Read data:
	std::vector<T> result(dims[0]);
	H5Dread(did, h5type<T>::get(), H5S_ALL, H5S_ALL, H5P_DEFAULT, result.data());
	H5Dclose(did);
	return result;
}

template<int n, typename T> std::vector< vectorn<n,T> > h5readMatrix(MPI_Comm comm, hid_t fid, const char* dname)
{	hid_t did = H5Dopen1(fid, dname);
	if(did<0){ fprintf(stderr, "Could not read dataset '%s' from input file.\n", dname); MPI_Abort(comm, 1); }
	//Get data dimensions:
	hid_t dspace = H5Dget_space(did);
	assert(H5Sget_simple_extent_ndims(dspace) == 2);
	hsize_t dims[2];
	H5Sget_simple_extent_dims(dspace, dims, NULL);
	assert(dims[1] == n);
	//Read data:
	std::vector< vectorn<n,T> > result(dims[0]);
	H5Dread(did, h5type<T>::get(), H5S_ALL, H5S_ALL, H5P_DEFAULT, result.data());
	H5Dclose(did);
	return result;
}

template<int m, int n, typename T> std::vector< matrixmn<m,n,T> > h5readMatrix(MPI_Comm comm, hid_t fid, const char* dname)
{	hid_t did = H5Dopen1(fid, dname);
	if(did<0){ fprintf(stderr, "Could not read dataset '%s' from input file.\n", dname); MPI_Abort(comm, 1); }
	//Get data dimensions:
	hid_t dspace = H5Dget_space(did);
	assert(H5Sget_simple_extent_ndims(dspace) == 3);
	hsize_t dims[3];
	H5Sget_simple_extent_dims(dspace, dims, NULL);
	assert(dims[1] == m);
	assert(dims[2] == n);
	//Read data:
	std::vector< matrixmn<m,n,T> > result(dims[0]);
	H5Dread(did, h5type<T>::get(), H5S_ALL, H5S_ALL, H5P_DEFAULT, result.data());
	H5Dclose(did);
	return result;
}

template<typename T> std::vector< std::vector<T> > h5readVectorArray(MPI_Comm comm, hid_t fid, const char* dname)
{	hid_t did = H5Dopen1(fid, dname);
	if(did<0){ fprintf(stderr, "Could not read dataset '%s' from input file.\n", dname); MPI_Abort(comm, 1); }
	//Get data dimensions:
	hid_t dspace = H5Dget_space(did);
	assert(H5Sget_simple_extent_ndims(dspace) == 2);
	hsize_t dims[2];
	H5Sget_simple_extent_dims(dspace, dims, NULL);
	//Divide outer dimension amongst processes:
	int nProcs, iProc;
	MPI_Comm_size(comm, &nProcs);
	MPI_Comm_rank(comm, &iProc);
	size_t oStart = (iProc * dims[0])/nProcs;
	size_t oStop = ((iProc+1) * dims[0])/nProcs;
	hsize_t count[] = { oStop-oStart, dims[1] };
	hsize_t offset[] = { oStart, 0 };
	//Initialize hyperslab:
	hid_t sid = H5Screate_simple(2, count, NULL);
	H5Sselect_hyperslab(dspace, H5S_SELECT_SET, offset, NULL, count, NULL);
	//Initialize collective access:
	hid_t plid = H5Pcreate(H5P_DATASET_XFER);
	H5Pset_dxpl_mpio(plid, H5FD_MPIO_COLLECTIVE);
	//Read data into temporary buffer:
	std::vector<T> buf(count[0]*count[1]);
	H5Dread(did, h5type<T>::get(), sid, dspace, plid, buf.data());
	//Copy from buffer to output structure:
	std::vector< std::vector<T> > result(dims[0]);
	for(size_t oDiff=0; oDiff<count[0]; oDiff++)
	{	size_t o = offset[0] + oDiff;
		result[o].resize(dims[1]);
		memcpy(result[o].data(), buf.data()+oDiff*dims[1], dims[1]*sizeof(T));
	}
	//Cleanup:
	H5Sclose(dspace);
	H5Sclose(sid);
	H5Dclose(did);
	H5Pclose(plid);
	return result;
}

template<int n, typename T> std::vector<std::vector<vectorn<n,T> > > h5readMatrixArray(MPI_Comm comm, hid_t fid, const char* dname)
{	hid_t did = H5Dopen1(fid, dname);
	if(did<0){ fprintf(stderr, "Could not read dataset '%s' from input file.\n", dname); MPI_Abort(comm, 1); }
	//Get data dimensions:
	hid_t dspace = H5Dget_space(did);
	assert(H5Sget_simple_extent_ndims(dspace) == 3);
	hsize_t dims[3];
	H5Sget_simple_extent_dims(dspace, dims, NULL);
	assert(dims[2]==n);
	//Divide outer dimension amongst processes:
	int nProcs, iProc;
	MPI_Comm_size(comm, &nProcs);
	MPI_Comm_rank(comm, &iProc);
	size_t oStart = (iProc * dims[0])/nProcs;
	size_t oStop = ((iProc+1) * dims[0])/nProcs;
	hsize_t count[] = { oStop-oStart, dims[1], dims[2] };
	hsize_t offset[] = { oStart, 0, 0 };
	//Initialize hyperslab:
	hid_t sid = H5Screate_simple(3, count, NULL);
	H5Sselect_hyperslab(dspace, H5S_SELECT_SET, offset, NULL, count, NULL);
	//Initialize collective access:
	hid_t plid = H5Pcreate(H5P_DATASET_XFER);
	H5Pset_dxpl_mpio(plid, H5FD_MPIO_COLLECTIVE);
	//Read data into temporary buffer:
	std::vector<vectorn<n,T> > buf(count[0]*count[1]);
	H5Dread(did, h5type<T>::get(), sid, dspace, plid, buf.data());
	//Copy from buffer to output structure:
	std::vector< std::vector<vectorn<n,T> > > result(dims[0]);
	for(size_t oDiff=0; oDiff<count[0]; oDiff++)
	{	size_t o = offset[0] + oDiff;
		result[o].resize(dims[1]);
		memcpy(result[o].data(), buf.data()+oDiff*dims[1], dims[1]*sizeof(vectorn<n,T>));
	}
	//Cleanup:
	H5Sclose(dspace);
	H5Sclose(sid);
	H5Dclose(did);
	H5Pclose(plid);
	return result;
}

template<typename T> void h5writeVector(MPI_Comm comm, hid_t fid, const char* dname, const std::vector<T>& data)
{	hid_t dataType = h5type<T>::get();
	//Create dataset:
	const hsize_t dims[] = { data.size() };
	hid_t sid = H5Screate_simple(1, dims, NULL);
	hid_t did = H5Dcreate(fid, dname, dataType, sid, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	H5Sclose(sid);
	if(did<0){ fprintf(stderr, "Could not create dataset '%s' in output file.\n", dname); MPI_Abort(comm, 1); }
	//Create collective write property:
	hid_t plid = H5Pcreate(H5P_DATASET_XFER);
	H5Pset_dxpl_mpio(plid, H5FD_MPIO_COLLECTIVE);
	//Write data:
	H5Dwrite(did, dataType, H5S_ALL, H5S_ALL, plid, data.data());
	//Cleanup:
	H5Dclose(did);
	H5Pclose(plid);
}

template<int n, typename T> void h5writeMatrix(MPI_Comm comm, hid_t fid, const char* dname, const std::vector< vectorn<n,T> >& data)
{	hid_t dataType = h5type<T>::get();
	//Create dataset:
	const hsize_t dims[] = { data.size(), n };
	hid_t sid = H5Screate_simple(2, dims, NULL);
	hid_t did = H5Dcreate(fid, dname, dataType, sid, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	H5Sclose(sid);
	if(did<0){ fprintf(stderr, "Could not create dataset '%s' in output file.\n", dname); MPI_Abort(comm, 1); }
	//Create collective write property:
	hid_t plid = H5Pcreate(H5P_DATASET_XFER);
	H5Pset_dxpl_mpio(plid, H5FD_MPIO_COLLECTIVE);
	//Write data:
	H5Dwrite(did, dataType, H5S_ALL, H5S_ALL, plid, data.data());
	//Cleanup:
	H5Dclose(did);
	H5Pclose(plid);
}

template<typename T> void h5getArrayDistribution(MPI_Comm comm, const char* dname, const std::vector<std::vector<T> >& data, size_t& oMin, size_t& oCount, size_t& oCountTot, size_t& iDim)
{	//Get inner and outer dimensions:
	oMin=data.size(); size_t oMax=0; //outer dimension extents (inclusive)
	iDim = 0; //inner dimension
	oCountTot = data.size();
	for(size_t o=0; o<oCountTot; o++)
		if(data[o].size())
		{	iDim = std::max(iDim, data[o].size());
			oMin = std::min(oMin, o);
			oMax = std::max(oMax, o);
		}
	//Check contiguity and uniformity:
	for(size_t o=oMin; o<=oMax; o++)
		if(iDim != data[o].size())
		{	int iProc; MPI_Comm_rank(comm, &iProc);
			fprintf(stderr, "Dataset '%s' is not contiguous and uniform on process %d.\n", dname, iProc);
			MPI_Abort(comm, 1);
		}
	unsigned long iDimMax = iDim;
	MPI_Allreduce(MPI_IN_PLACE, &iDimMax,1, MPI_UNSIGNED_LONG, MPI_MAX, comm);
	if(!iDim) iDim = iDimMax;
	unsigned long iDimMin = iDim;
	MPI_Allreduce(MPI_IN_PLACE, &iDimMin,1, MPI_UNSIGNED_LONG, MPI_MIN, comm);
	if(iDimMin != iDimMax)
	{	fprintf(stderr, "Dataset '%s' inner dimensions not constant across processes (range: %lu-%lu).\n", dname, iDimMin, iDimMax);
		MPI_Abort(comm, 1);
	}
	//Check total size:
	oCount = oMax>=oMin ? oMax-oMin+1 : 0;
	unsigned long oCountSum = oCount; MPI_Allreduce(MPI_IN_PLACE, &oCountSum,1, MPI_UNSIGNED_LONG, MPI_SUM, comm);
	if(oCountTot != oCountSum)
	{	fprintf(stderr, "Dataset '%s' per-process dimension sum (%lu) does not match total (%lu).\n", dname, oCountSum, oCountTot);
		MPI_Abort(comm, 1);
	}
}

template<typename T> void h5writeVectorArray(MPI_Comm comm, hid_t fid, const char* dname, const std::vector<std::vector<T> >& data)
{	hid_t dataType = h5type<T>::get();
	//Determine total and per-process data dimensions:
	size_t oMin, oCount, oCountTot, iDim;
	h5getArrayDistribution(comm, dname, data, oMin, oCount, oCountTot, iDim);
	const hsize_t dims[] = { oCountTot, iDim }; //total dataset dimensions
	const hsize_t count[] = { oCount, iDim }; //local dataset dimensions
	const hsize_t offset[] = { oMin, 0 }; //offset of local dataset into total
	//Copy data into consecutive memory space:
	std::vector<T> dataFlat(count[0]*count[1]);
	for(size_t oDiff=0; oDiff<count[0]; oDiff++)
		memcpy(dataFlat.data()+oDiff*count[1], data[oMin+oDiff].data(), count[1]*sizeof(T));
	//Create dataset:
	hid_t sidGlobal = H5Screate_simple(2, dims, NULL);
	hid_t did = H5Dcreate(fid, dname, dataType, sidGlobal, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	H5Sclose(sidGlobal);
	if(did<0){ fprintf(stderr, "Could not create dataset '%s' in output file.\n", dname); MPI_Abort(comm, 1); }
	//Select hyperslab:
	hid_t sid = H5Screate_simple(2, count, NULL);
	sidGlobal = H5Dget_space(did);
	H5Sselect_hyperslab(sidGlobal, H5S_SELECT_SET, offset, NULL, count, NULL);
	//Create collective write property:
	hid_t plid = H5Pcreate(H5P_DATASET_XFER);
	H5Pset_dxpl_mpio(plid, H5FD_MPIO_COLLECTIVE);
	//Write data:
	if(!oCount)
	{	H5Sselect_none(sid);
		H5Sselect_none(sidGlobal);
	}
	H5Dwrite(did, dataType, sid, sidGlobal, plid, dataFlat.data());
	//Cleanup:
	H5Sclose(sidGlobal);
	H5Sclose(sid);
	H5Dclose(did);
	H5Pclose(plid);
}

template<int n, typename T> void h5writeMatrixArray(MPI_Comm comm, hid_t fid, const char* dname, const std::vector<std::vector<vectorn<n,T> > >& data)
{	hid_t dataType = h5type<T>::get();
	//Determine total and per-process data dimensions:
	size_t oMin, oCount, oCountTot, iDim;
	h5getArrayDistribution(comm, dname, data, oMin, oCount, oCountTot, iDim);
	const hsize_t dims[] = { oCountTot, iDim, n }; //total dataset dimensions
	const hsize_t count[] = { oCount, iDim, n }; //local dataset dimensions
	const hsize_t offset[] = { oMin, 0, 0 }; //offset of local dataset into total
	//Copy data into consecutive memory space:
	std::vector< vectorn<n,T> > dataFlat(count[0]*count[1]);
	for(size_t oDiff=0; oDiff<count[0]; oDiff++)
		memcpy(dataFlat.data()+oDiff*count[1], data[oMin+oDiff].data(), count[1]*sizeof(vectorn<n,T>));
	//Create dataset:
	hid_t sidGlobal = H5Screate_simple(3, dims, NULL);
	hid_t did = H5Dcreate(fid, dname, dataType, sidGlobal, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	H5Sclose(sidGlobal);
	if(did<0){ fprintf(stderr, "Could not create dataset '%s' in output file.\n", dname); MPI_Abort(comm, 1); }
	//Select hyperslab:
	hid_t sid = H5Screate_simple(3, count, NULL);
	sidGlobal = H5Dget_space(did);
	H5Sselect_hyperslab(sidGlobal, H5S_SELECT_SET, offset, NULL, count, NULL);
	//Create collective write property:
	hid_t plid = H5Pcreate(H5P_DATASET_XFER);
	H5Pset_dxpl_mpio(plid, H5FD_MPIO_COLLECTIVE);
	//Write data:
	if(!oCount)
	{	H5Sselect_none(sid);
		H5Sselect_none(sidGlobal);
	}
	H5Dwrite(did, dataType, sid, sidGlobal, plid, dataFlat.data());
	//Cleanup:
	H5Sclose(sidGlobal);
	H5Sclose(sid);
	H5Dclose(did);
	H5Pclose(plid);
}

#endif //TETRAHEDRONTRANSPORT_H5IO_H
