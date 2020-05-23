classdef RasterVariables<handle
    methods(Access = public,Sealed=true)
        function var=Initialize(obj)
            % the function must be called after the basinMask has been
            % obtained
            [rows,columns]=size(obj.basinMask);
            var=zeros(rows,columns);
            var(~obj.basinMask)=NaN;
        end
        function [ind,bInBasin]=sub2indInBasin(obj,matRow,matCol)
            % ind: one dimensional indices that are within the basin
            % bInBasin: logical vector that indicates whether (matRow,matCol) is within the basin
            bInBasin=InRectangle(obj,matRow,matCol);
            ind=zeros(length(matRow),1);
            ind(bInBasin)=sub2ind(size(obj.basinMask),matRow(bInBasin),matCol(bInBasin));
%             indexOut=obj.basinMask(ind);
%             ind(~indexOut)=[];
            ind(bInBasin)=ind(bInBasin).*obj.basinMask(ind(bInBasin));
            bInBasin(bInBasin)=ind(bInBasin)>0;
            ind(~bInBasin)=NaN;
        end
        function index=InRectangle(obj,matRow,matCol)
            [rows,columns]=size(obj.DEM);
            index=logical((matRow<=rows).*(matRow>0).*(matCol<=columns).*(matCol>0));
        end
        
        function SaveRaster(obj,mat,fileGeoTif)
            geotiffwrite(fileGeoTif,mat,obj.geo);
        end
    end
    methods(Static=true,Access=protected)
        function [value,bDistributed]=readVarInfo(gfileID,keywordType,keywordVar,commentSymbol)
            value=-1;
            bDistributed=-1;
            if ~isempty(keywordType)
                while bDistributed==-1
                    tline = fgetl(gfileID);
                    if strcmp(tline(1),commentSymbol)==1
                        continue;
                    end
                    strArr = regexp(tline,commentSymbol,'split');
                    strContent=strArr{1};
                    strContent=strtrim(strContent);
                    if ~isempty(strfind(strContent, keywordType))
                        strValue=regexp(strContent,'=','split');
                        if strcmpi(strValue{2},'distributed')==1
                           bDistributed=true;
                        else
                           bDistributed=false;
                        end
                    end
                end
            else
                bDistributed=true;
            end
            while value==-1
                tline = fgetl(gfileID);
                if strcmp(tline(1),commentSymbol)==1
                    continue;
                end
                strArr = regexp(tline,commentSymbol,'split');
                strContent=strArr{1};
                strContent=strtrim(strContent);
                if ~isempty(strfind(strContent, keywordVar))
                    strValue=regexp(strContent,'=','split');
                    if ~isnan(str2double(strValue{2}))
                       value=str2double(strValue{2});
                    else
                       value=strValue{2};
                    end
                end 
            end
        end
        function bRes=IsGCS(refMat)
            if(refMat(2,1)<2)
                bRes=true;
            else
                bRes=false;
            end
        end
    end
    properties(Access=public)
        basinMask;%the logical index matrix of valid data elements
        spatialRef;%OsGeo.OGR.SpatialReference
        geoTrans; % geographic transformation coefficients. Used to convert(r,c) to (mapX,mapY) or(lat lon)
        bGCS % indicate whether the projection of the basic data is gcs or PCS
        proj;% projection info if has
    end
    methods(Abstract,Static)
        fileNames = GenerateFileNames(obj,dirFolder)
    end
end