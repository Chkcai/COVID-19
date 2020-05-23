function PrepareObsAndGlobal()
streamRas='C:\Data\Basin data\poyang Lake\Ganjiang\export\stream.asc';
basinRas='C:\Data\Basin data\poyang Lake\Ganjiang\export\mask.asc';
obsDir='C:\Data\Basin data\poyang Lake\observations';
site_shp={'C:\Data\Basin data\poyang Lake\HydroSites\sites_river_off_correct.shp';...
          'C:\Data\Basin data\poyang Lake\HydroSites\sites_in_river.shp'};
StreamflowData='C:\Data\Basin data\poyang Lake\gauge_data.xlsx';
BasinDir='C:\simulation\Ganjiang';
gName='Ganjiang_Multisites.Project';
NoData=-9999;
outlet=62302250;
[site_no,lat,lon,nRow,nCol,bottom,left,res]=GetSitesInBasin(site_shp,streamRas,basinRas);
years=[1997,1998,1999,2000,2001,2002,2003,2004,2005];
sNULL=ExtractRunoff(StreamflowData,years,site_no,obsDir,NoData);
site_no(sNULL)=[];
lat(sNULL)=[];
lon(sNULL)=[];
GenerateGlobalCTR(BasinDir,gName,bottom,left,res,nRow,nCol,NoData,site_no,lat,lon,outlet);
end