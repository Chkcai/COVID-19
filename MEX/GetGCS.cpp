//calling convention
//wkt=GetGCS(strCode);
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
    char* pszShortName=NULL;
    char* pszWKT=NULL;
    mxArray *mxProjection;
    OGRSpatialReference oSRS;
    if ( nrhs !=1)
		mexErrMsgTxt("incorrect input"); 
    pszShortName=ImportString((mxArray *)prhs[0]);
    oSRS.SetWellKnownGeogCS(pszShortName);
    oSRS.exportToWkt( &pszWKT );
    mxProjection=mxCreateString(pszWKT);
    plhs[0]=mxProjection;
}