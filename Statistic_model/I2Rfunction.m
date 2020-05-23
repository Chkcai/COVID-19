function r=I2Rfunction(data,g)
% Estimate reproductive number from infection number
dbstop if error
[m,~]=size(data);
r=NaN(m,1);
for i=2:m
    ct=data(i);
    s=0;
    s2=0;
    if i>20
        for j=i-20:i-1
                s=s+data(j)*g(i-j);
                s2=s2+g(i-j);
        end
    else
        for j=1:i-1
                s=s+data(j)*g(i-j);
                s2=s2+g(i-j);
        end
    end
    ct_1=s/s2;
    r(i)=ct/ct_1;
end
