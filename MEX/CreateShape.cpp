/* GDAL MEX Utility (v1.0) by Shen,Xinyi
 *updated on Mar, 2020 to allow other types than points
   contact: Xinyi.Shen@uconn.edu,Jan, 2016
   read a shape file */
/*calling convention
*CreateShape(shapeName,type,spatialRef,fieldName1,fieldType1, fieldName2,fieldType2,...)
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
    /*******************input********************/
    const char *shapeName, *fieldName;
    int geoType;
    OGRwkbGeometryType geoTypeGDAL;
    OGRFieldType fieldType;
    const char *spatialRef;
    const char *pszDriverName = "ESRI Shapefile";
    OGRSFDriverH hDriver;
    shapeName=ImportString(prhs[0]);
    geoType=(int)*mxGetPr(prhs[1]);
    spatialRef=ImportString(prhs[2]);
    
    // CREATE DATASOURCE
    hDriver = OGRGetDriverByName( pszDriverName );
    if( hDriver == NULL )
    {
        printf( "%s driver not available.\n", pszDriverName );
        exit( 1 );
    }
    OGRDataSourceH hDS;
    OGRSpatialReferenceH srH=OSRNewSpatialReference (spatialRef);
    hDS = OGR_Dr_CreateDataSource( hDriver, shapeName, NULL );
    
    if( hDS == NULL )
    {
        printf( "Creation of output file failed.\n" );
        exit( 1 );
    }
    // create Layer
    OGRLayerH hLayer;
    switch (geoType)
    {
        case 1:
            geoTypeGDAL=wkbPoint;
            break;
        case 2:
            geoTypeGDAL=wkbPolygon;
            break;
    }
    
    hLayer = OGR_DS_CreateLayer( hDS, "", srH, geoTypeGDAL, NULL );
    if( hLayer == NULL )
    {
        printf( "Layer creation failed.\n" );
        exit( 1 );
    }
    
    // create Fields
    int nFields=(nrhs-3)/2;
    if (nFields<=0)
        printf( "Empty shape file!\n");
    for(int ifield=1;ifield<=nFields;ifield++)
    {
        fieldName=ImportString(prhs[2*ifield+1]);
        fieldType=(OGRFieldType)(int)(*mxGetPr(prhs[2*ifield+2]));
        OGRFieldDefnH hFieldDefn;
        hFieldDefn = OGR_Fld_Create( fieldName, fieldType );
        //OGR_Fld_SetWidth( hFieldDefn, 32);
        if( OGR_L_CreateField( hLayer, hFieldDefn, TRUE ) != OGRERR_NONE )
        {
            printf( "Creating %s field failed.\n",fieldName );
            exit( 1 );
        }
        OGR_Fld_Destroy(hFieldDefn);
    }
//  
    OGR_DS_Destroy( hDS );
}