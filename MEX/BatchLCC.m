function BatchLCC(dirIn,dirOut,projNew)
curFile = mfilename('fullpath');
[MEXDir,~,~]=fileparts(curFile);
[MDir,~,~]=fileparts(MEXDir);
addpath([MDir,'/','blockwise']);
files = dir(dirIn);
noDataSrc=[0,256];
noDataDst=0;
blockSize=6000;
for i=1:length(files)
    fileIn=files(i).name;
    if strcmpi(fileIn,'..') || strcmpi(fileIn,'.')
        continue;
    end
    disp(['processing ' fileIn num2str(i) ' of ' num2str(length(files))])
    ifileIn=[dirIn,fileIn];
    fileTemp=[dirOut,'temp.tif'];
    fileOut=[dirOut,fileIn];
    %% compute the target resolution
    [~,~,~,geoTrans,proj]=RasterInfo(ifileIn);
    [Y1,X1]=RowCol2Proj(geoTrans,[1;1;2],[1;2;1]);
    [X2,Y2]=ProjTransform(proj,projNew,X1,Y1);
    resXNew=X2(2)-X2(1);
    resYNew=Y2(3)-Y2(1);
    resNew=(resXNew-resYNew)/2;
    resXNew=resNew;
    resYNew=-resNew;
    imProjTrans(ifileIn,projNew,resXNew,resYNew,fileTemp);
    SetNodata(fileTemp,fileOut,noDataSrc,noDataDst,blockSize);
    delete(fileTemp);
end
end