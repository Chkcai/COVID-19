function [area,gridArea]=RasterArea(raster,geoTrans,bGCS)
[rows,columns]=size(raster);
basinMask=~isnan(raster);
r=1:rows;
c=1:columns;
[C,R]=meshgrid(c,r);
R=R(basinMask);
C=C(basinMask);
%calculate the area of each grid
if bGCS
    LenSN=abs(geoTrans(6))*110574.0;
    [lat,~]=RowCol2Proj(geoTrans,R,C);
    LenEW=LenSN*cosd(lat);
else
    LenSN=abs(geoTrans(6));
    LenEW=abs(geoTrans(2));
end
A=LenSN*LenEW*1e-6; % Convert to km^2
gridArea=nan(rows,columns);
gridArea(basinMask)=A;
area=sum(gridArea(basinMask));
end