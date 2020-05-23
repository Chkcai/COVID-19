function res=C2Ifunction(data)
dbstop if error
for i=1:14
    a(i,1)=logncdf(i,1.621,0.418);
    % A lognormal distribution is employed by with model;
end
w(1)=a(1);
for i=2:14
    w(i,1)=a(i)-a(i-1);
end
[m,~]=size(data);
s=zeros(m,1);
for i=1:m
    if data(i)<0
        data(i)=0;
    end
    if i>=14
        in_w=flipud(w);
        nIn=data(i)*in_w/sum(w);
        s(i-13:i,:)=s(i-13:i,:)+nIn;
    else
        in_w_s=w(1:i);
        in_w=flipud(in_w_s);
        nIn=data(i)*in_w;
        s(1:i,:)=s(1:i,:)+nIn/sum(w);
    end
end
res=s;
end