function copyVRT(fileList,dirLoc,rasterName)
fid=fopen(fileList,'r');
fileListLoc=[dirLoc,'fileList.txt'];
fidLoc=fopen(fileListLoc,'w');
fileName=[];
while fileName~=-1
    fileName=fgetl(fid);
    fprintf(fidLoc,[fileName,'\n']);
    [~,name,ext]=fileparts(fileName);
    switch ext
        case '.img'
            copyHFA(fileName,[dirLoc,name,ext]);
        otherwise
            copyfile(fileName,[dirLoc,name,ext]);
    end
end
fclose(fidLoc);
fclose(fid);
cmd=[' !gdalbuildvrt -input_file_list ',fileListLoc, ' ', dirLoc,rasterName];
eval(cmd);
end