/* GDAL MEX Utility (v1.0) by Shen,Xinyi
   contact: Xinyi.Shen@uconn.edu,Feb, 2015
   coordinates transformation*/ 
/************************Calling Convention************************/
/* [X2,Y2]=ProjTransform(strWKT1,strWKT2,X1,Y1)           */
/******************************************************************/
#include "mexOperation.h"
#include "gdal.h"
#include "ogr_srs_api.h"
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
    /*************input variables********************************/
    char *strWkt1, *strWkt2;
    double *X1, *Y1;
    /************ output variables******************************/
    double *X2, *Y2;
    /***********************************************************/
    int nCount,M,N,size;
    OGRCoordinateTransformationH coordianteTransformationH;
    OGRSpatialReferenceH srH1,srH2;
    /*************Converting input variables********************/
    if (nrhs!=4)
        mexErrMsgTxt("incorrect output");
    strWkt1=ImportString(prhs[0]);
    strWkt2=ImportString(prhs[1]);
    X1=mxGetPr(prhs[2]);
    Y1=mxGetPr(prhs[3]);
    M=mxGetM(prhs[2]);
    N=mxGetN(prhs[2]);
    nCount=M*N;
    X2=new double[nCount];
    Y2=new double[nCount];
    memcpy(X2,X1,nCount*sizeof(double));
    memcpy(Y2,Y1,nCount*sizeof(double));
    /**********************************************************/
    srH1=OSRNewSpatialReference (strWkt1);
    srH2=OSRNewSpatialReference (strWkt2);
    coordianteTransformationH=OCTNewCoordinateTransformation(srH1,srH2);
    OCTTransform(coordianteTransformationH,nCount,X2,Y2,NULL);
    OCTDestroyCoordinateTransformation(coordianteTransformationH);
    size = GDALGetDataTypeSize ( GDT_Float64 ) / 8;
    plhs[0]=ExportRealMatrix(X2,M, N,size,false);
    plhs[1]=ExportRealMatrix(Y2,M, N,size,false);
    
}