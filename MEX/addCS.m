function addCS(xllcorner,yllcorner,cellsize,nRows,dataType,NoDataVal,fileUnProj,fileProjInfo)
if isempty(fileProjInfo)
    fileProjInfo=fileUnProj;
end
[~,~,~,~,proj]=RasterInfo(fileProjInfo);
if ~isempty(xllcorner) && ~isempty(yllcorner) && ~isempty(cellsize) && isempty(nRows)
    geoTrans=zeros(1,6);
    geoTrans(1)=xllcorner+cellsize/2;
    geoTrans(2)=cellsize;
    geoTrans(3)=0;
    geoTrans(4)=yllcorner+(nRows-1/2)*cellsize;
    geoTrans(5)=0;
    geoTrans(6)=-cellsize;
else
    [~,~,~,geoTrans]=RasterInfo(fileProjInfo);
end
if isempty(dataType)
    [~,~,~,~,~,dataType]=RasterInfo(fileProjInfo);
end
if isempty(NoDataVal)
    [~,~,~,~,~,~,NoDataVal]=RasterInfo(fileProjInfo);
end
outFormat='GTiff';
WriteRasterInfo(fileUnProj,geoTrans,proj,dataType,outFormat,NoDataVal);
% WriteRaster(fileProj,raster,geoTrans,proj,dataType,outFormat,NoDataVal)
end