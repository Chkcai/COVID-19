function trim_err_rec(fileIn,filecsvOut,dtStart,dtEnd,region,colsIndep,varargin)
%% This function is used to preprocess the data in the regression mode.
% Triming days or regions with poor data quality
% fmt='%{d/mm/yyyy} %d %f %f %f %f %f %f %f %f %d %d %f %f';
useDeath=false;
opt='trans_ratio';
minCase=30;
tMinTrans=20;
NPrior=0;N0Prior=0;
relative_time=true;
inputWindow=0;
for i = 1 : 2 : length(varargin)
    switch varargin{i}
        case 'comm_prior'
            comm_prior=varargin{i+1};
            NPrior=comm_prior(1);
            N0Prior=comm_prior(2);
        case 'relative_time'
            relative_time=varargin{i+1};
        case 'useDeath'
            useDeath=varargin{i+1};
        case 'optY'
            opt=varargin{i+1};
        case 'minCase'
            minCase=varargin{i+1};
        case 'tMinTrans'
            tMinTrans=varargin{i+1};
        case 'input'
            inputArg=varargin{i+1};
            inputWindow=inputArg(1);
            incubation=inputArg(2);
    end
end
T=readtable(fileIn,'sheet',3);
objCode=T{:,'OBJ'};
dates=T{:,'Date'};
confirm=T{:,'Confirmed'};
dConfirm=T{:,'delta'};
if useDeath
    R=T{:,'ratio_DC'};
else
    R=T{:,'ratio'};
end
if iscell(R)
    R=char(R);
    R=str2num(R);
end
dn=datenum(dates);
indepMat=T{:,colsIndep};
% mask=isnan(mat(:,1)) ;
mask=dn<dtStart | dn>dtEnd;
dn(mask)=[];
objCode(mask)=[];
confirm(mask)=[];
R(mask)=[];
indepMat(mask,:)=[];
dConfirm(mask)=[];
%% get the number of dates and regions
udn=unique(dn);
nDates=length(udn);
nRegions=length(unique(objCode));
%% locate the minimum day for community spreading and take the value as N0
minDay=accumarray(objCode,confirm,[],@(x)findMinDay(x,minCase),[],true);
N=reshape(confirm,[nDates,nRegions]);
objCode=reshape(objCode,[nDates,nRegions]);
dnMat=reshape(dn,[nDates,nRegions]);
dnMatAbs=dnMat;
dConfirm=reshape(dConfirm,[nDates,nRegions]);
[uCode,~,minDay]=find(minDay);
col=uCode-uCode(1)+1;
noCom=isinf(minDay);
minDay(noCom)=nDates;


switch opt
    case 'growth_exp'
        ind=sub2ind([nDates,nRegions],minDay,col);
        N0=N(ind)';
        N0(noCom)=Inf;
        nDays=(1:nDates)';
        matNDays=repmat(nDays,[1,nRegions]);
        dDays=matNDays-repmat(minDay',[nDates,1]);
        Y=(log(N)-log(repmat(N0,[nDates,1])))./dDays;
        %% remove the days before community spreading
        for iC=1:nRegions
            Y(1:minDay(iC),iC)=NaN;
            dnMat(minDay(iC)+1:end,iC)=dnMat(minDay(iC)+1:end,iC)-dnMat(minDay(iC),iC);
            dnMat(1:minDay(iC),iC)=NaN;
        end
    case 'growth'
        Y=(N(2:end,:)-N(1:end-1,:))./N(1:end-1,:);
        index=false(nDates,nRegions);
        index(1,:)=true;
        index=index(:);
        for iC=1:nRegions
            Y(1:minDay(iC)-1,iC)=NaN;
            dnMat(minDay(iC)+1:end,iC)=dnMat(minDay(iC)+1:end,iC)-dnMat(minDay(iC),iC);
            dnMat(1:minDay(iC),iC)=NaN;
        end
        objCode(1,:)=[];
        dnMat(1,:)=[];
        dConfirm(1,:)=[];
        indepMat(index,:)=[];
    case 'trans_ratio'
        Y=reshape(R,[nDates,nRegions]);
        if relative_time
            for iC=1:nRegions
                Y(1:minDay(iC),iC)=NaN;
                dnMat(minDay(iC)+1:end,iC)=dnMat(minDay(iC)+1:end,iC)-dnMat(minDay(iC),iC);
                dnMat(1:minDay(iC),iC)=NaN;
            end
        else
            dnMat=dnMat-min(min(dnMat))+1;
        end
end
%% trim off non-community spreading days if the N-prior days has less than N0 nonzero newly confirmed days
% offset dConfirm to remove zeros before the community spread time
for iC=1:nRegions
    dConfirm(1:minDay(iC)-1,iC)=max(1,dConfirm(1:minDay(iC)-1,iC));
end
if N0Prior>0 && NPrior>N0Prior
    hasZeroNew=dConfirm<=0;
    template=ones(NPrior,1);
    for iC=1:nRegions
        zeroNum=conv(template,hasZeroNew(:,iC));
        zeroNum(nDates+1:end)=[];
        indEnd=find(zeroNum>N0Prior,1,'first');
        iIndEnd=find(indEnd>minDay(iC),1,'first');
        indEnd=indEnd(iIndEnd);
        if ~isempty(indEnd)
            Y(indEnd:end,iC)=NaN;
        end
    end
end
if inputWindow>0
    hWin=round(inputWindow/2);
    index0=1:nDates;
    indexCompare0=index0-hWin;
    indexCompare0(indexCompare0<1)=1;
    indexCompare0((indexCompare0+inputWindow)>nDates)=nDates-inputWindow+1;
    ratio=nan(nDates,inputWindow);
    for iC=1:nRegions
        for iW=1:inputWindow
            ratio(:,iW)=dConfirm(:,iC)./dConfirm(indexCompare0+iW-1,iC);
        end
        sortedRatio=sort(ratio,2,'ascend');
        minRatio=sortedRatio(:,2);
        indInput=find(minRatio>10);
        for iInput=1:length(indInput)
            Y(max(indInput(iInput)-incubation+1,1):indInput(iInput),iC)=NaN;
        end
    end
end
%% trim off regions with too short community spreading days
nEff=sum(~isnan(dnMat) & (~isnan(Y)),1);
col2Remove=nEff<tMinTrans;
Y(:,col2Remove)=NaN;
Y=Y(:);
dnMat=dnMat(:);
dnMatAbs=datestr(dnMatAbs(:),'mm/dd/yyyy');
objCode=objCode(:);
regionNo=region*ones(length(dnMat),1);
mat=[dnMat,Y,regionNo,double(objCode),indepMat];
% mat(isnan(Y),:)=[];
C=[cellstr(dnMatAbs),num2cell(mat)];
clear mat
% dates=datestr(dn);
% C=[cellstr(dates),C];
T = cell2table(C,'VariableNames',{'Date_Abs','Date',opt,'region','OBJ',colsIndep{:}});
% T = cell2table(C,'VariableNames',{'Date',opt,'OBJ'...
%     'Precipitation','Relative_humidity','Specific_humidity',...
%     'Solar_radiation','UV_radiation','NO2','Temperature','wind'});
writetable(T,filecsvOut);
end
function ind=findMinDay(cases,tMinCase)
ind=find(cases>tMinCase,1,'first');
if isempty(ind)
    ind=Inf;                                                                                                                                                                                                                                        
end
end