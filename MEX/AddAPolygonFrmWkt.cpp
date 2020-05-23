/* GDAL MEX Utility (v1.0) by Shen,Xinyi
   contact: Xinyi.Shen@uconn.edu,Jan, 2016
   read a shape file */
/*calling convention
*AddAPolygon(shapeName,wktPolygon,fieldName1,fieldValue1, fieldName2,fieldValue2, ...)
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
    #pragma comment (lib,"gdal_i.lib")sa
    #define strEqual 0
    #define strcmp strcmp
#else
    #error Platform not supported
    #define strEqual 0
#endif
void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] )
{
    const char *shapeName, *fieldName, *fieldValue;
    char* wktString;
    const char *pszDriverName = "ESRI Shapefile";
    int nPoints;
    double* xLoc, *yLoc;
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
    wktString=ImportString(prhs[1]);
    // Create a Feature from user input
    OGRFeatureH hFeature;
    OGRGeometryH hPolygon;
    // create a feature from definition
    hFeature = OGR_F_Create( OGR_L_GetLayerDefn( hLayer ) );
    //Fill fields
    int nFields=(nrhs-2)/2;
    for(int ifield=1;ifield<=nFields;ifield++)
    {
        fieldName=ImportString(prhs[2*ifield]);
        fieldValue=ImportString(prhs[2*ifield+1]);
        OGR_F_SetFieldString( hFeature, OGR_F_GetFieldIndex(hFeature, fieldName), fieldValue );
    }
    
    OGR_G_CreateFromWkt(&wktString,NULL,&hPolygon);
    OGR_F_SetGeometry( hFeature, hPolygon ); 
    OGR_G_DestroyGeometry(hPolygon);

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