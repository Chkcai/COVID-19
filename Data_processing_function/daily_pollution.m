function daily_pollution(dirInRas,shapeAdmin,dirList,shapeOut,dtStart,dtEnd,fileGrid,varargin)
%% This function computes the daily average of a certain pollution concentration
% a repensentative (mean) value of each polygon will be output
%% input
% dirInRas: the directory of All L2 product with daily subfolders in yyyy-mm-dd format
% shapeAdmin: polygons with administrative boundary
% dirList: the directory of the searched result overpassing this directory
% shapeOut: the output name
% dtStart,dtEnd: start and ending time of the computation in datenumber
% fileGrid: the regular grid raster file whoes spatial resolution is close to the product
%%
qac=true;
for i = 1 : 2 : length(varargin)
    switch varargin{i}
        case 'qa'
            qac=varargin{i+1};
    end
end
%% main
NoData=9.969209968386869e+36; % from gdalinfo
thQA=50; % https://sentinel.esa.int/documents/247904/2474726/Sentinel-5P-Level-2-Product-User-Manual-Nitrogen-Dioxide Section 8.6 
OS=computer;
if strcmpi(OS,'GLNXA64')
    pathSplitor='/';
elseif strcmpi(OS,'PCWIN64')
    pathSplitor='\';
end

curFile = mfilename('fullpath');
[curDir,~,~]=fileparts(curFile);
[MfileDir,~,~]=fileparts(curDir);
addpath([MfileDir,'/','MEX']);
GDALLoad();

% fileGrid=[dirTemp,'grid.tif'];
%% read the administrative regions. If the output already exists, read it instead and append new dates
shapeRegion=shapeAdmin;
if exist(shapeOut,'file')==2
    shapeRegion=shapeOut;
end
S=shaperead(shapeRegion);
% adminCode=[S.OBJECTID];
adminCode=[S.OBJECT_ID];
% adminCode=str2num(char(adminCode));
% %% if the regular grids do not exist, generate it for once
% if exist(fileGrid,'file')~=2
%    cmd=['! gdal_rasterize -a ST_CNTY_CO -burn 1 -of GTiff -tr 0.04 0.04 ',...
%             shapeRegion, ' ',fileGrid];
%    eval(cmd);
% end
[gridIndex,geoTrans]=ReadRaster(fileGrid);
[rows,cols]=size(gridIndex);
maskGrid=~isnan(gridIndex);
[rowGrid,colGrid]=find(gridIndex);
indAllGrids=sub2ind([rows,cols],rowGrid(maskGrid),colGrid(maskGrid));
nGrids=max(indAllGrids);
%% loop over all dates to ouput the administrative mean
fmt='%s %s';
for dtCur=dtStart:dtEnd
    dtStr=datestr(dtCur,'yyyy-mm-dd');
    dirDay=[dirInRas,dtStr,pathSplitor];
    if ~isempty(dirList)
        fileList=dir([dirList,'*',dtStr,'*.csv']);
        fileList=[fileList.folder,pathSplitor,fileList.name];
        fid=fopen(fileList,'r');
        C=textscan(fid,fmt,'Delimiter',',');
        swath=C{1};
        fclose(fid);
    end
%     swath=dir([dirDay,swath,'*.nc']);
    NO2_grid=sparse(nGrids,1);
    count=sparse(nGrids,1);
    varName=['NO2_',datestr(dtCur,'yymmdd')];
    fprintf(1,'processing on %s\n',dtStr)
    for iSwath = 1:length(swath)
        % read the nc file
        fprintf(1,'processing on swatch %d\n',iSwath)
%         fileS5P=dir([dirDay,'*',swath(iSwath).name]);
        fileS5P=dir([dirDay,'*',swath{iSwath},'*.nc']);
        fileS5P=[dirDay,fileS5P.name];
        if exist(fileS5P,'file')~=2
            disp(['Do not have ', fileS5P, '. Skip...']);
            continue;
        end
        fileLat=['HDF5:"',fileS5P,'"://PRODUCT/latitude'];
        fileLon=['HDF5:"',fileS5P,'"://PRODUCT/longitude'];
        fileQA=['HDF5:"',fileS5P,'"://PRODUCT/qa_value'];
        fileNO2=['HDF5:"',fileS5P,'"://PRODUCT/SUPPORT_DATA/DETAILED_RESULTS/nitrogendioxide_summed_total_column'];
        %% prevent a corrupted netcdf to crash the program
        cmdTestNCCorruption=['! gdalinfo ',fileS5P];
        fileInfo=evalc(cmdTestNCCorruption);
        if contains(fileInfo,'Segmentation fault','IgnoreCase',true)
            disp('NetCDF file corrputed. Skip...')
            continue;
        end
        %% end of the prevention
        try
            NO2=ReadRaster(fileNO2);
            lat=ReadRaster(fileLat);
            lon=ReadRaster(fileLon);
            if qac
                qa=ReadRaster(fileQA);
                % remove Nodata
                mask=NO2==NoData  | NO2<0 | qa<thQA ;
            else
                mask=NO2==NoData  | NO2<0;
            end
            NO2(mask)=[];lat(mask)=[];lon(mask)=[];
            % regularize all grids
            [row,col]=Proj2RowCol(geoTrans,lat,lon,false);
            mask=row<1 | row>rows | col<1 | col>cols;
            NO2(mask)=[];row(mask)=[];col(mask)=[];
            % average over each grid
            ind=sub2ind([rows,cols],row,col);
            % remove grids out of the region
            objID=gridIndex(ind);
            mask=isnan(objID);
            ind(mask)=[];
            NO2(mask)=[];
            iNO2_grid=accumarray(ind',NO2',[],@mean,[],true);
            iCount=accumarray(ind',ones(length(ind),1),[],@sum,[],true);
            [ind,~,iNO2_grid]=find(iNO2_grid);
            [~,~,iCount]=find(iCount);
            NO2_grid(ind)=(iNO2_grid.*iCount+NO2_grid(ind).*count(ind))./(iCount+count(ind));
            count(ind)=count(ind)+iCount;
        catch
            continue;
        end
    end
    
    [ind,~,NO2_grid]=find(NO2_grid);
    objID=gridIndex(ind);
    NO2_Admin=accumarray(objID,NO2_grid,[],@mean,[],true);
    [objID,~,NO2_Admin]=find(NO2_Admin);
    [uAdminCode,~,ic]=unique(adminCode);
    [isContained,LoB]=ismember(objID,uAdminCode);
    NO2_Admin(~isContained)=[];
    LoB(~isContained)=[];
    uVal=nan(length(uAdminCode),1);
    uVal(LoB)=NO2_Admin;
    uVal(isnan(uVal))=-1;% Some region has no data
    val=num2cell(uVal(ic));% for some administrative file, a unique id is repeated
    [S.(varName)]=val{:};
%     LoA=ismember(adminCode,objID);
%     NoValue=-1*ones(sum(~LoA),1);
%     NoValue=num2cell(NoValue);
%     [S(~LoA).(varName)]=NoValue{:};
    shapewrite(S,shapeOut);
end

end