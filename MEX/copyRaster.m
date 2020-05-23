function copyRaster(fileSrc,fileDst,fmtDst,addNodata,blockSize,showProg)
%% This function copies a raster file by reading and writing it
[nR,nC,r0,c0,nRows,nCols]=GetBlock(fileSrc,blockSize);
[~,RowImage,ColImage,geoTrans,proj,dataType,NodataVal]=RasterInfo(fileSrc);
for iR=1:nR
    for iC=1:nC
        if showProg
            prog=((iR-1)*nC+iC)/(nR*nC);
            fprintf(1,'%d%%...',int32(prog*100));
        end
        rasIn=ReadRaster(fileSrc,r0(iR),c0(iC),nRows(iR),nCols(iC));
        if ~isempty(addNodata)
            rasIn(rasIn==addNodata)=NaN;
        end
        WriteRaster(fileDst,rasIn,geoTrans,proj,dataType,fmtDst,NodataVal,r0(iR),c0(iC),RowImage,ColImage);
    end
end
end