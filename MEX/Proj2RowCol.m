function [row,col]= Proj2RowCol(geoTrans, projY, projX,keepFrac)
%% GDAL MEX Utility (v1.0) by Shen,Xinyi
%   contact: Xinyi.Shen@uconn.edu,Feb, 2015
%   (Ygeo,Xgeo)->(row,col)
% this function is the inverse of RowCol2Proj and can be verified by solving (row,col) from (geoX,geoY)
dTemp = geoTrans(2)*geoTrans(6) - geoTrans(3)*geoTrans(6);  
col = (geoTrans(6)*(projX - geoTrans(1)) - geoTrans(3)*(projY - geoTrans(4))) / dTemp + 0.5;  
row = (geoTrans(2)*(projY - geoTrans(4)) - geoTrans(5)*(projX - geoTrans(1))) / dTemp + 0.5;   
% (row,col) system in matlab is (1,1) based
if ~exist('keepFrac','var') || ~keepFrac
    col = round(col);  
    row = round(row);
end
end