/* GDAL MEX Utility (v1.0) by Shen,Xinyi
   contact: Xinyi.Shen@uconn.edu,Feb, 2015
   read raster information from an image file */

//calling convention
//[nBands,Rows,Cols,geoTrans,proj,dataType,NodataVal]=RasterInfo(fileName);
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
    /* input variables*/
	mxArray* mxSrcFile;
	int strLen,status;
	char* pszSrcFile;
	GDALDatasetH pDSrc;
    GDALRasterBandH pBand;
    int CPLErr,bSuccess;
    int xSize,ySize,nBands;
    //geoTransform
    double geoTrans[6];mxArray *mxGeoTrans;
    double NoData;mxArray* mxNoData;
    //spatial Reference
    const char *spatialRef;mxArray *mxProjection;
    // data type
    GDALDataType DataType;mxArray *mxDataType;
    
    int out_type_size;
    void* dptr;
    mxArray* mxInput[2];
	if ( nrhs !=1)
		mexErrMsgTxt("incorrect input"); 
    pszSrcFile=ImportString((mxArray *)prhs[0]);
    //Start GDAL Operations
	pDSrc = (GDALDatasetH)GDALOpen(pszSrcFile, GA_ReadOnly);
    mxFree(pszSrcFile);
    pBand=GDALGetRasterBand(pDSrc,1);
    DataType=GDALGetRasterDataType(pBand);
    xSize=GDALGetRasterXSize(pDSrc);
    ySize=GDALGetRasterYSize(pDSrc);
    nBands=GDALGetRasterCount(pDSrc);
   
    if (nlhs<1)
        mexErrMsgTxt ( "output arguments not assigned" );
    else
    {
        plhs[0] = mxCreateDoubleScalar(nBands);
        if (nlhs>=3)
        {
            plhs[1]=mxCreateDoubleScalar(ySize);
            plhs[2]=mxCreateDoubleScalar(xSize);
        }
        if (nlhs>=4)
        {
            CPLErr=GDALGetGeoTransform(pDSrc,geoTrans);
            if (CPLErr!=CE_None)
                mexErrMsgTxt("error in reading head"); 
            out_type_size = GDALGetDataTypeSize ( GDT_Float64 ) / 8;
            mxGeoTrans=ExportRealMatrix((const double*)geoTrans,1,6,out_type_size,0);
            plhs[3]=mxGeoTrans;
        }
        if (nlhs>=5)
        {
            spatialRef=GDALGetProjectionRef (pDSrc);
            mxProjection=mxCreateString(spatialRef);
            plhs[4]=mxProjection;
        }
        if (nlhs>=6)
        {
            mxDataType=mxCreateDoubleScalar(DataType);
            plhs[5]=mxDataType;
        }
        if (nlhs>=7)
        {
            NoData=GDALGetRasterNoDataValue(pBand,&bSuccess);
            mxNoData=mxCreateDoubleScalar(NoData);
            plhs[6]=mxNoData;
        }
    }
    GDALClose((GDALDatasetH) pDSrc); 
}