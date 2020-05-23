function ERA5_read(Input_path,Output_path)
%% This function generate tif from ERA5 grib file.

% Input
% Input_path: the directory of all ERA5 grib file.
% Output_path: the directory of the generated tif file 
%% Load
addpath([MfileDir,'/','MEX']);
GDALLoad()
if strcmpi(OS,'GLNXA64')
    pathSplitor='/';
elseif strcmpi(OS,'PCWIN64')
    pathSplitor='\';
end
Input_path=[Input_path pathSplitor];
Output_path=[Output_path pathSplitor];
%% Main
grb_file=dir([Input_path '*.grib']);
for i=1:length(grb_file)
    file_name=grb_file(i,:).name;
    [nBands,~,~,geoTrans,proj,dataType,NodataVal]=RasterInfo([Input_path file_name]);
    nDate=nBands/24;
    Out_file=[Output_path file_name(1:end-5)];
    if exist(Out_file,'dir')~=7
        mkdir(Out_file);
    end
    for i_date=1:nDate
        for i_hour=1:24
            this_band=i_hour+24*(i_date-1);
            [raster,~,~,~,~]=ReadMultiBandRaster([Input_path file_name],this_band);
            if i_date<10
                str_date=['0' num2str(i_date)];
            else
                str_date= num2str(i_date);
            end
            if i_hour<10
                str_hour=['0' num2str(i_hour)];
            else
                str_hour=num2str(i_hour);
            end
            fileRas=[Out_file pathSplitor str_date str_hour '.tif'];
            WriteRaster(fileRas,raster,geoTrans,proj,dataType,'GTiff',NodataVal(this_band))
            disp([file_name(1:end-5) str_date str_hour ' finished']);
            [row,col]=size(raster);
            Ind_ra=sub2ind(size(raster),1:row*col);
            data_hourly(:,i_hour)=raster(Ind_ra);
        end
        daily_res=NaN(row,col);
        daily_res(Ind_ra)=nanmean(data_hourly(Ind_ra,:),2);
        fileRas_d=[Out_file pathSplitor str_date '_mean.tif'];
        WriteRaster(fileRas_d,daily_res,geoTrans,proj,dataType,'GTiff',NodataVal(this_band))
        clear data_hourly
    end
end
end