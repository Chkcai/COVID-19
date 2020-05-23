function copyHFA(fileSrc,fileDst)
[~,~,extSrc]=fileparts(fileSrc);
[~,~,extDst]=fileparts(fileDst);
if strcmpi(extSrc,'.img') && strcmpi(extDst,'.img')
    fileSrcIGE=fileSrc;
    fileSrcIGE(end-3:end)='.ige';
    fileDstIGE=fileDst;
    fileDstIGE(end-3:end)='.ige';
    copyfile(fileSrc,fileDst);
    if exist(fileSrcIGE,'file')==2
        copyfile(fileSrcIGE,fileDstIGE);
    end
else
    disp('not ERDAS image format. No files are copied')
end

end