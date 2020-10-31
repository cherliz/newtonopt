function [df] = logisticDF(x0)
a=x0(1); b=x0(2); c=x0(3);

[time,data]=runAll();
time = zeros(length(data),1);
for i=1:1:length(data)
    time(i) = i;
end

[row,col] = size(data);
df = zeros(1,3);
err = zeros(3,1);
erra = zeros(3,1);
errb = zeros(3,1);
errc = zeros(3,1);

for i=1:1:row
    err(i) = data(i)-(a/(1+b*exp(-1*c+time(i))));
    erra = -1/(1+b*exp(-1*c*time(i)));
    errb = a*exp(-1*c*time(i))/((1+b*exp(-1*c*time(i)))^2);
    errc = -1*a*b*time(i)*exp(-1*c*time(i))/((1+b*exp(-1*c*time(i)))^2);
end

for i = 1:1:row
    df(1) = df(1) + 2*err(i)*erra(i);
    df(2) = df(2) + 2*err(i)*errb(i);
    df(3) = df(3) + 2*err(i)*errc(i); 
end

end