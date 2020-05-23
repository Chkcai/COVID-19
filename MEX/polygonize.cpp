#include "mexOperation.h"
// calling convention
// polygonize(strSrcRas,strDstFile)
// #ifndef  WIN32
// #include <unistd.h>
// #endif
#include "gdal.h"
//#include "cpl_conv.h"
#include "gdal_priv.h"
#include "gdal_alg.h"
#include "mex.h"
#include "ogr_api.h"
#include "ogrsf_frmts.h"
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

void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] )
{
    const char *pszSrcFile, *pszDstFile; 
    const char *pszFormat="ESRI Shapefile";  
    //fileName
    pszSrcFile=ImportString(prhs[0]);
    pszDstFile=ImportString(prhs[1]);
    GDALDataset* poSrcDS = (GDALDataset*) GDALOpen(pszSrcFile, GA_ReadOnly);
    OGRSFDriver* poDriver = (OGRSFDriver*)OGRSFDriverRegistrar::GetRegistrar()->GetDriverByName(pszFormat);
    if (poDriver == NULL)
    {
        GDALClose((GDALDatasetH)poSrcDS);
        mexErrMsgTxt("Cannot create the output file type.");
        return;
    }
    OGRDataSource* poDstDS = poDriver->CreateDataSource(pszDstFile, NULL);
    if (poDstDS == NULL)
    {
       GDALClose((GDALDatasetH)poSrcDS);
       mexErrMsgTxt("Cannot create the output file");
       return;
    }
    OGRSpatialReference *poSpatialRef = new OGRSpatialReference(poSrcDS->GetProjectionRef());
//    const char *strLayerName = OGRLayer::GetName(pszDstFile); 
    OGRLayer* poLayer = poDstDS->CreateLayer("raster", poSpatialRef, wkbPolygon, NULL);
    if (poLayer == NULL)
    {
        GDALClose((GDALDatasetH)poSrcDS);
        OGRDataSource::DestroyDataSource(poDstDS);
        delete poSpatialRef;
        poSpatialRef = NULL;
        mexErrMsgTxt("Cannot create layer");
        return;
    }
    OGRFieldDefn ofieldDef("ID", OFTInteger);    //create field table
    if (poLayer->CreateField(&ofieldDef) != OGRERR_NONE)
    {
        GDALClose((GDALDatasetH)poSrcDS);
        OGRDataSource::DestroyDataSource(poDstDS);
        delete poSpatialRef;
        poSpatialRef = NULL;
        mexErrMsgTxt("failed to create field table");
        return;
    }
    GDALRasterBandH hSrcBand = (GDALRasterBandH) poSrcDS->GetRasterBand(1); 
    if (GDALPolygonize(hSrcBand, NULL, (OGRLayerH)poLayer, 0, NULL, NULL, NULL) != CE_None)
    {
        GDALClose((GDALDatasetH)poSrcDS);
        OGRDataSource::DestroyDataSource(poDstDS);
        delete poSpatialRef;
        poSpatialRef = NULL;
        mexErrMsgTxt("failed");
        return;
    }
    GDALClose((GDALDatasetH)poSrcDS); 
    OGRDataSource::DestroyDataSource(poDstDS);  
    delete poSpatialRef;
    poSpatialRef = NULL;
}