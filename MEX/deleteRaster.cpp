/* GDAL MEX Utility (v1.0) by Shen,Xinyi
   contact: Xinyi.Shen@uconn.edu,Feb, 2015
   remove a raster file */
// deleteRaster(fileName,format)
#include "mexOperation.h"
#include "gdal.h"
#include "mex.h"
#ifdef __linux__
    #pragma comment (lib,"libgdal.so")
    //linux code goes here
#elif _WIN64
    // windows code goes here
    #pragma comment (lib,"gdal_i.lib")
#else
    #error Platform not supported
#endif


void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] )
{
    const char* format, *fileName;
    GDALDriverH driverH;
    fileName=ImportString(prhs[0]);
    format=ImportString(prhs[1]);
    driverH = GDALGetDriverByName(format);
    GDALDeleteDataset(driverH,fileName);
}