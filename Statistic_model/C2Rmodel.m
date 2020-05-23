function [this_data_delta,this_Infection,this_ratio]=C2Rmodel(data)
%% 
%This model estimates the daily reproductive number from daily cumulative 
%confirmed number 

% Input
% data:Daily cumulative confirmed number data(in a row/col);

% Output
%this_data_delta: Daily confirmed number;
%this_Infection: Daily infected number estimated by the statistic model,
%i.e.equation (1);
%this_ratio: Daily reproductive number;
%% Main
nDate=length(data);
%The CDF of infection(I) to reproductive number(R)
load('serial_interval.mat');
g=d(:,2);%Load from the mat
[sRow,~]=find(~isnan(data)&data~=0,1);% Find the date of first confirmed
if ~isempty(sRow)
    this_data=data(sRow:end,1);
    this_data_delta(1)=this_data(1);
    for j=2:length(this_data)
        this_data_delta(j,1)=this_data(j)-this_data(j-1);
        % Delta of confirmed number
    end
    this_Infection=C2Ifunction(this_data_delta);%equation (1)
    this_ratio=I2Rfunction(this_Infection,g);%equation (2)
    r_lag=NaN(sRow-1,1);% Add NaN to the date without confirmed number
    this_Infection=[r_lag;this_Infection];
    this_ratio=[r_lag;this_ratio];
    this_data_delta=[r_lag;this_data_delta];
else% If no confirmed number in this region
    this_Infection=NaN(size(data));
    this_ratio=NaN(size(data));
    this_data_delta=NaN(size(data));
end
end
