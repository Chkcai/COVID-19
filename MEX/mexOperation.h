#include "mex.h"
#include "gdal.h"
#include "ogr_spatialref.h"
#include <math.h>
#include <cstring>
#include <string.h>
#pragma comment (lib,"gdal_i.lib")

#ifdef __linux__
    #pragma comment (lib,"libgdal.so")
    #define strEqual 0
    #define strcmp strcmp
    //linux code goes here
#elif _WIN64
    // windows code goes here
    #pragma comment (lib,"gdal_i.lib")
    #define strEqual 1
    #define strcmp strcmp
#else
    #error Platform not supported
    #define strEqual 0
#endif
//Transpose a matrix
mxArray * Transpose(mxArray * mxMat1)
{
    const mwSize *dims;
    mxArray *mxMat2;
    dims=mxGetDimensions(mxMat1);
    mxMat2 = mxCreateNumericArray(2, dims, mxDOUBLE_CLASS, mxREAL );
    mexCallMATLAB(1,&mxMat2,1, &mxMat1,"transpose");
    return mxMat2;
}
//Copy the matlab string to C++ string
char* ImportString(const mxArray * mxStr)
{
    char *pszStr;
    int strLen,status;
	if ( mxIsChar(mxStr) != 1 ) 
	{
		mexErrMsgTxt("Input file name must be a string\n" );
	}
	strLen = mxGetN(mxStr) + 1; 
	pszStr = (char *) mxCalloc ( strLen, sizeof(char) );
	status = mxGetString ( mxStr, pszStr, strLen );
	if ( status != 0 ) 
	{
		mexErrMsgTxt ( "sting copy error.\n" );
	}
    return pszStr;
}
mxArray* ExportRealMatrix(const double* raster,int rows,int cols,int size,bool bTranspose)
{
    mwSize *dims;
    void * dptr;
    mxArray *mxMat1;
    mxArray *mxMat2;
    dims=new mwSize;
    if (bTranspose)
    {
        dims[0]=cols;
        dims[1]=rows;    
    }
    else
    {
        dims[1]=cols;
        dims[0]=rows;
    }
    mxMat1 = mxCreateNumericArray(2, (const mwSize*)dims, mxDOUBLE_CLASS, mxREAL );
    dptr = mxGetPr ( mxMat1 );
    memcpy (dptr, raster, size*cols*rows );
    if (bTranspose)
    {
        mxMat2=Transpose(mxMat1);
        mxDestroyArray(mxMat1);
    }
    else
        mxMat2=mxMat1;
    return mxMat2;
}

bool exists(const char* fileName)
{
    mxArray *mxFileName;
    mxArray *indicator;
    mxArray *mxInput[2];
    int result;
    mxFileName=mxCreateString(fileName);
    indicator=mxCreateDoubleScalar(0);
    mxInput[0]=mxFileName;
    mxInput[1]=mxCreateString("file");
    mexCallMATLAB(1,&indicator,1,mxInput,"exist");
    result=(bool)(*mxGetPr(indicator));
    return result;
}
char* fileExtension(const char* fullPath)
{
    mxArray *mxDir,*mxFileName,*mxExt,*mxFullPath;
    mxArray *fileOut[3];
    mxFullPath=mxCreateString(fullPath);
    fileOut[0]=mxDir;
    fileOut[1]=mxFileName;
    fileOut[2]=mxExt;
    mexCallMATLAB(3,fileOut,1,&mxFullPath,"fileparts");
    return ImportString(fileOut[2]);
}
bool IsGeographic(const char *wktSrc, double *geoTrans)
{
    bool bGeo;
    OGRSpatialReferenceH srSrcH;
    srSrcH=OSRNewSpatialReference (wktSrc);
    if (OSRIsGeographic(srSrcH))
        bGeo=true;
    else if (OSRIsProjected(srSrcH))
        bGeo=false;
    else// if the information of spatialReference cannot tell the type of projection
        // the resolution will be used for judgement
        // if the resolution is larger than 5(m), it's a projected
        // this rule requires modification
    {
        if (fabs(geoTrans[1])>5)
            bGeo=false;
        else
            bGeo=true;
    }
    return bGeo;
}