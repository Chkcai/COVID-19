function PopCal(popInRas,NO2InShp,fileGrid,Temp_path,shapeOut,UpLevelInShp)
%% This function estimates the population of each administrative boundary 
% used in NO2 result, and gives the attribution in the higher 
% administrative boundary 

% Input
% popInRas: the directory of population density data in tif format
% NO2InShp: the directory of the NO2 shp
% fileGrid: the regular grid raster file whoes spatial resolution is same
% to the population data.
% Temp_path: the directory of temporary files
% shapeOut: 
% UpLevelInShp: polygons of higher level administrative boundary
%% Main
OS=computer;
if strcmpi(OS,'GLNXA64')
    pathSplitor='/';
elseif strcmpi(OS,'PCWIN64')
    pathSplitor='\';
end
addpath '/shared/stormcenter/Shen/src/MEX-2.3.0';
GDALLoad();
%% Read Object_id
S=shaperead(NO2InShp);
adminCode=[S.OBJECT_ID];
[gridIndex,geoTrans,proj]=ReadRaster(fileGrid);
[rows,cols]=size(gridIndex);
[rowGrid,colGrid]=find(~isnan(gridIndex));
indAllGrids=sub2ind([rows,cols],rowGrid,colGrid);
%% Resample the Population
if exist(Temp_path,'dir')~=7
    mkdir(Temp_path);
end
temp_pop_tif=[Temp_path pathSplitor 'Temp_pop.tif'];
ResampleAndClip(geoTrans,proj,cols,rows,popInRas,temp_pop_tif,'GTiff',1,1);
Pop_val=ReadRaster(temp_pop_tif);
Ojb_id=gridIndex(indAllGrids);
Pop=Pop_val(indAllGrids);
[mask,~]=find(Pop==0|isnan(Pop));
Pop_grid=accumarray(Ojb_id,Pop);
[uAdminCode,~,ic]=unique(adminCode);
[objID,~]=find(Pop_grid~=0);
[~,LoB]=ismember(objID,uAdminCode);
uVal=nan(length(uAdminCode),1);
uVal(LoB)=Pop_grid;
uVal(isnan(uVal))=-1;% Some region has no data
val=num2cell(uVal(ic));% for some administrative file, a unique id is repeated
varName='POP';
[S.(varName)]=val{:};
%% Upscale to new region
province_l2={S.NAME_1};
S_tar=shaperead(UpLevelInShp);
province_l1={S_tar.NAME_1};
OBJ=NaN(length(province_l2),1);
for i=1:length(province_l1)
    province_l1_str=char(province_l1(i));
    str_len=length(province_l1_str);
    for j=1:length(province_l2)
        province_l2_str=char(province_l2(j));
        if strcmp(province_l1_str,province_l2_str)
            OBJ(j)=i;
        end
    end
end
varName_U='UpObj';
val_U=num2cell(OBJ);
[S.(varName_U)]=val_U{:};
shapewrite(S,shapeOut);
