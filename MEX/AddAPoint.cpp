/* GDAL MEX Utility (v1.0) by Shen,Xinyi
   contact: Xinyi.Shen@uconn.edu,Jan, 2016
   read a shape file */
/*calling convention
*AddAPoint(shapeName,xLoc,yLoc,fieldName1,fieldValue1, fieldName2,fieldValue2, ...)
*/
#include "mexOperation.h"
#include "gdal.h"
#include "ogr_api.h"
#ifdef __linux__
    #pragma comment (lib,"libgdal.so")
    #define strEqual 0
    #define strcmp strcmp
    //linux code goes here
#elif _WIN64
    // windows code goes here
    #pragma comment (lib,"gdal_i.lib")
    #define strEqual 0
    #define strcmp strcmp
#else
    #error Platform not supported
    #define strEqual 0
#endif
void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] )
{
    const char *shapeName, *fieldName, *fieldValue;
    const char *pszDriverName = "ESRI Shapefile";
    double xLoc,yLoc;
    shapeName=ImportString(prhs[0]);
    // Open a existing shapefile to write
    OGRLayerH hLayer;
    OGRSFDriverH hDriver;
    OGRDataSourceH hDS;
    hDriver = OGRGetDriverByName( pszDriverName );
    if( hDriver == NULL )
    {
        printf( "%s driver not available.\n", pszDriverName );
        exit( 1 );
    }
    
    hDS = OGROpen( shapeName, TRUE, NULL );
    if( hDS == NULL )
    {
        printf( "Open failed.\n" );
        exit( 1 );
    }
    hLayer = OGR_DS_GetLayer( hDS, 0 );
    if( hLayer == NULL )
    {
        printf( "Layer creation failed.\n" );
        exit( 1 );
    }
    // read location from user input
    xLoc=(double)(*mxGetPr(prhs[1]));
    yLoc=(double)(*mxGetPr(prhs[2]));
    // Create a Feature from user input
    OGRFeatureH hFeature;
    OGRGeometryH hPt;
    // create a feature from definition
    hFeature = OGR_F_Create( OGR_L_GetLayerDefn( hLayer ) );
    //Fill fields
    int nFields=(nrhs-3)/2;
    for(int ifield=1;ifield<=nFields;ifield++)
    {
        fieldName=ImportString(prhs[2*ifield+1]);
        fieldValue=ImportString(prhs[2*(ifield+1)]);
        OGR_F_SetFieldString( hFeature, OGR_F_GetFieldIndex(hFeature, fieldName), fieldValue );
    }
    hPt = OGR_G_CreateGeometry(wkbPoint);
    OGR_G_SetPoint_2D(hPt, 0, xLoc, yLoc);
    
    OGR_F_SetGeometry( hFeature, hPt ); 
    OGR_G_DestroyGeometry(hPt);

    if( OGR_L_CreateFeature( hLayer, hFeature ) != OGRERR_NONE )
    {
       printf( "Failed to create the feature in shapefile.\n" );
       exit( 1 );
    }
    // dispose the feature
    OGR_F_Destroy( hFeature );
    // dispose the data source
    OGR_DS_Destroy( hDS );
}