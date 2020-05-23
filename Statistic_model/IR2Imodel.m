function this_Infection=IR2Imodel(data)
%%
% This model estimates the infection number from daily infection number
% with a new daily reproductive number.
% This function was used to estimate the Global transmissibility reduction 
% through restricting SA in the paper.

% Input: Data is a matrix with m row and 2 col. The first col is the daily
% infection number, and the second col is the new daily reproductive
% number.

% Output
% this_Infection: Daily infection number under the new daily reproductive
% number.
%% Main
[nDate,~]=size(data);
load('serial_interval.mat');%d
%The CDF of infection(I) to reproductive number(R)
g=d(:,2);% Load from the mat
[sRow,~]=find(~isnan(data(:,1))&data(:,1)~=0,1);
if ~isempty(sRow)
    this_data=data(sRow:end,:);
    this_new_inf=R2Ifunction(this_data(:,1),g,this_data(:,2));
    r_lag=NaN(sRow-1,1);
    this_Infection=[r_lag;this_new_inf];
else
    this_Infection=NaN(size(data));
end
end