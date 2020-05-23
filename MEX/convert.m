function convert()
GDALLoad();
curFile = mfilename('fullpath');
[curDir,~,~]=fileparts(curFile);
[progDir,~,~]=fileparts(curDir);
addpath([progDir,'/blockwise']);
fileIn='/shared/manoslab/data/GEEBAM_v2p2_dNBR_Classes.tif';
fileOut='/shared/manoslab/data/GEEBAM_v2p2_dNBR_Classes_unique.tif';
[nR,nC,r0,c0,nRows,nCols,rowsImage,colsImage]=GetBlock(fileIn,3000,0);
dataType=1;
fmtDst='GTiff';
NodataVal=255;
for iR=1:nR
    for iC=1:nC
        prog=((iR-1)*nC+iC)/(nR*nC);
        fprintf(1,'%d%%...',int32(prog*100));
        [rasIn,geoTrans,proj]=ReadRaster(fileIn,r0(iR),c0(iC),nRows(iR),nCols(iC));
        rasOut=nan(nRows(iR),nCols(iC));
        rasOut(rasIn>=8 & rasIn<=9)=4;
        rasOut(rasIn>=6 & rasIn<=7)=3;
        rasOut(rasIn>=4 & rasIn<=5)=2;
        rasOut(rasIn>=1 & rasIn<=3)=1;
        WriteRaster(fileOut,rasOut,geoTrans,proj,dataType,fmtDst,NodataVal,r0(iR),c0(iC),rowsImage,colsImage);
    end
end
end