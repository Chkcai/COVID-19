function crop(fileSrc,fileRegion,fileDst)
%% crop the region defined by fileRegion from fileSrc. 
% the projection of fileSrc and fileRegion can be different
% fileDst is the cropped file. 
% Different than Resample and clip, the cropping result, fileDst is of the
% same projection and resolution as fileSrc, instead of fileRegion
%% update history
% 1) Dec. 1 2019, Force the grids must be integer values of fileSrc
GDALLoad();
[~,rows,cols,geoTransRegion,projRegion]=RasterInfo(fileRegion);
[~,~,~,geoTransSrc,projSrc]=RasterInfo(fileSrc);
resXSrc=geoTransSrc(2);
resYSrc=geoTransSrc(6);
[YCorner,XCorner]=RowCol2Proj(geoTransRegion,[1,rows],[1,cols]);
[XNewCorner,YNewCorner]=ProjTransform(projRegion,projSrc,[XCorner(1),XCorner(1),XCorner(2),XCorner(2)]',[YCorner(1),YCorner(2),YCorner(1),YCorner(2)]');
XLT=min(XNewCorner);
YLT=max(YNewCorner);
XRB=max(XNewCorner);
YRB=min(YNewCorner);

%% begin 1)
[rowLT,colLT]=Proj2RowCol(geoTransSrc,YLT,XLT,true);
colLT=floor(colLT);
rowLT=ceil(rowLT);
[YLT,XLT]=RowCol2Proj(geoTransSrc,rowLT,colLT);
%% end 1)
geoTransNew=zeros(6,1);
geoTransNew(1)=XLT-resXSrc/2;
geoTransNew(4)=YLT-resYSrc/2;
geoTransNew(2)=resXSrc;
geoTransNew(6)=resYSrc;
%% begin 1)
rowsNew=ceil((YRB-YLT)/resYSrc)+1;
colsNew=ceil((XRB-XLT)/resXSrc)+1;
ResampleAndClip(geoTransNew,projSrc,colsNew,rowsNew,fileSrc,fileDst,'GTiff',1,1);
% rowsNew=(YRB-YLT)/resYSrc+1;
% colsNew=(XRB-XLT)/resXSrc+1;
% ResampleAndClip(geoTransNew,projSrc,colsNew,rowsNew,fileSrc,fileDst,'GTiff',1);
%% end 1)
end