function daily_meteorologyVariable(dirInRas,popInRas,shapeAdmin,fileGrid,shapeOut)
%% This function computes the daily average of a meteorology variables with the weight of population
% a repensentative (mean) value of each polygon will be output
%% input
% dirInRas: the directory of daily meteorology variable in tif format
% popInRas: the directory of population density data in tif format
% shapeAdmin: polygons with administrative boundary
% shapeOut: 
% fileGrid: the regular grid raster file whoes spatial resolution is same
% to the population data.
%% Main
varname='Var';
OS=computer;
if strcmpi(OS,'GLNXA64')
    pathSplitor='/';
elseif strcmpi(OS,'PCWIN64')
    pathSplitor='\';
end
addpath([MfileDir,'/','MEX']);
GDALLoad();
%% read the administrative regions. If the output already exists, read it instead and append new dates
shapeRegion=shapeAdmin;
if exist(shapeOut,'file')==2
    shapeRegion=shapeOut;
end
S=shaperead(shapeRegion);
adminCode=[S.OBJECT_ID];
[gridIndex,geoTrans,proj]=ReadRaster(fileGrid);
[rows,cols]=size(gridIndex);
[rowGrid,colGrid]=find(~isnan(gridIndex));
indAllGrids=sub2ind([rows,cols],rowGrid,colGrid);
%% Population data resample
[outpath,~,~]=fileparts(shapeOut);
Temp_path=[outpath pathSplitor 'Temp'];
if exist(Temp_path,'dir')~=7
    mkdir(Temp_path);
end
temp_pop_tif=[Temp_path pathSplitor 'Temp_pop.tif'];
ResampleAndClip(geoTrans,proj,cols,rows,popInRas,temp_pop_tif,'GTiff',1,1);
%% loop over all dates
tif_file=dir([dirInRas '*_mean.tif']);
for i=1:length(tif_file)
    this_tif=[dirInRas tif_file(i,:).name];
    if i<10
        varName=[varname '0' num2str(i)];
    else
        varName=[varname num2str(i)];
    end
    [Var_tif,geoTrans_tif,proj_tif,dataType,NodataVal]=ReadRaster(this_tif);
    Var_trans=[Var_tif,Var_tif,Var_tif];
    trans_val=[360 0 0 0 0 0];
    geoTrans_new=geoTrans_tif-trans_val;
    Temp_trans_tif=[Temp_path pathSplitor 'temp_trans' num2str(i) '.tif'];
    WriteRaster(Temp_trans_tif,Var_trans,geoTrans_new,proj_tif,dataType,'GTiff',NodataVal)
    Temp_clip_tif=[Temp_path pathSplitor 'temp_clip' num2str(i) '.tif'];
    ResampleAndClip(geoTrans,proj,cols,rows,Temp_trans_tif,Temp_clip_tif,'GTiff',1,2);
    Var_val=ReadRaster(Temp_clip_tif);
    Pop_val=ReadRaster(temp_pop_tif);
    Ojb_id=gridIndex(indAllGrids);
    Var=Var_val(indAllGrids);
    Pop=Pop_val(indAllGrids);
    [mask,~]=find(isnan(Var)|Pop==0|isnan(Pop));
    Var(mask)=[];Pop(mask)=[];Ojb_id(mask)=[];
    Var_pop=Var.*Pop;
    Var_grid=accumarray(Ojb_id,Var_pop);
    Pop_grid=accumarray(Ojb_id,Pop);
    [uAdminCode,~,ic]=unique(adminCode);
    [objID,~]=find(Var_grid~=0&Pop_grid~=0);
    Mean_grid=Var_grid(objID)./Pop_grid(objID);
    [~,LoB]=ismember(objID,uAdminCode);
    uVal=nan(length(uAdminCode),1);
    uVal(LoB)=Mean_grid;
    uVal(isnan(uVal))=-1;% Some region has no data
    val=num2cell(uVal(ic));% for some administrative file, a unique id is repeated
    [S.(varName)]=val{:};
    shapewrite(S,shapeOut);
    disp([num2str(i) ' finished'])
end
rmdir(Temp_path,'s');
