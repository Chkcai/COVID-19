function cp_shape_head(shpDst,shpSrc)
[dirIn,fileIn]=fileparts(shpSrc);
[dirOut,fileOut]=fileparts(shpDst);
copyfile([dirIn,'\',fileIn,'.prj'],[dirOut,'\',fileOut,'.prj']);
end