function polygon=parsePolygon(strValue)
    strPos=strsplit(strValue,{'(((','((',',',')))','))'});
    if length(strPos)>1 && contains(strPos{1},'POLYGON','IgnoreCase',true)
        p1=strtrim(strPos{2});p1s=strsplit(p1,' ');
        p2=strtrim(strPos{3});p2s=strsplit(p2,' ');
        p3=strtrim(strPos{4});p3s=strsplit(p3,' ');
        p4=strtrim(strPos{5});p4s=strsplit(p4,' ');
        polygon=polyshape([str2double(p1s{1}) str2double(p2s{1})...
            str2double(p3s{1}) str2double(p4s{1})],...
            [str2double(p1s{2}) str2double(p2s{2})...
            str2double(p3s{2}) str2double(p4s{2})]);
    else
        polygon=[];
    end
end