function imProjTrans(fileImage,projNew,resXNew,resYNew,fileNew)
% this function entirely changes the CS of a image file
%% input
% fileImage : the orignal file with complete CS info
% projNew: projNew the target projection
% fileNew: the output file
GDALLoad();
[~,rows,cols,geoTrans,proj]=RasterInfo(fileImage);
[YCorner,XCorner]=RowCol2Proj(geoTrans,[1,rows],[1,cols]);
[XNewCorner,YNewCorner]=ProjTransform(proj,projNew,[XCorner(1),XCorner(1),XCorner(2),XCorner(2)]',[YCorner(1),YCorner(2),YCorner(1),YCorner(2)]');
XLT=min(XNewCorner);
YLT=max(YNewCorner);
XRB=max(XNewCorner);
YRB=min(YNewCorner);
geoTransNew=zeros(6,1);
geoTransNew(1)=XLT-resXNew/2;
geoTransNew(4)=YLT-resYNew/2;
geoTransNew(2)=resXNew;
geoTransNew(6)=resYNew;
rowsNew=(YRB-YLT)/resYNew+1;
colsNew=(XRB-XLT)/resXNew+1;
ResampleAndClip(geoTransNew,projNew,colsNew,rowsNew,fileImage,fileNew,'GTiff',1);
end