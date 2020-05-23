#include "mexOperation.h"
#include "gdal.h"
#include "ogr_spatialref.h"
#ifdef __linux__
    #pragma comment (lib,"libgdal.so")
    //linux code goes here
#elif _WIN64
    // windows code goes here
    #pragma comment (lib,"gdal_i.lib")
#else
    #error Platform not supported
#endif
//IsGeographic(spatialRef,geoTrans)
void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] )
{
    const char* wkt;
    double* geoTrans;
    bool bGeo;
    if (nrhs<2)
        mexErrMsgTxt("insufficient input");
    wkt=ImportString(prhs[0]);
    geoTrans=mxGetPr(prhs[1]);
    bGeo=IsGeographic(wkt, geoTrans);
    plhs[0]=mxCreateDoubleScalar(bGeo);
}


