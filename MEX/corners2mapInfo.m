function [XLT,YLT,cols,rows]=corners2mapInfo(proj1,proj2,resX,resY,left,right,top,btm)
[XLT,YLT]=ProjTransform(proj1,proj2,left,top);
[XRB,YRB]=ProjTransform(proj1,proj2,right,btm);
rows=round(abs((YRB-YLT)/resY))+1;
cols=round(abs((XRB-XLT)/resX))+1;
end