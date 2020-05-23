function Ninfection=R2Ifunction(data,g,r)
[m,~]=size(data);
Ninfection=NaN(m,1);
Ninfection(1)=data(1);
for i=2:m
    if ~isnan(r(i))
        s=0;
        s2=0;
        if i>20
            for j=i-20:i-1
                s=s+Ninfection(j)*g(i-j);
                s2=s2+g(i-j);
            end
        else
            for j=1:i-1
                s=s+Ninfection(j)*g(i-j);
                s2=s2+g(i-j);
            end
        end
        ct_1=s/s2;
        Ninfection(i)=ct_1*r(i);
    else
        Ninfection(i)=data(i);
    end
end
end