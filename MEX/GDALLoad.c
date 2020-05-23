#include "gdal.h"
#include "mex.h"
#ifdef __linux__
    #pragma comment (lib,"libgdal.so")
    /* linux code goes here */
#elif _WIN64
    /* windows code goes here */
    #pragma comment (lib,"gdal_i.lib")
#else
    #error Platform not supported
#endif
void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray*prhs[] )
{
	if ( nrhs > 0 || nlhs>0 )
	{
		mexErrMsgTxt("no argument is allowed."); 
	}
	GDALAllRegister();
	OGRRegisterAll();
}