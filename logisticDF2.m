function [df2] = logisticDF2(x0)
% this calculates Df2, so second derivative of least squares
% err and partial derivatives!
[x,data]=runAll();
time = zeros(length(data),1);
for i=1:1:length(data)
    time(i) = i;
end
a=x0(1); b=x0(2); c=x0(3);

a=x0(1); b=x0(2); c=x0(3);
[row,col] = size(data);
df2 = zeros(3);
err = zeros(3,1);
erra = zeros(3,1);
errb = zeros(3,1);
erraa = zeros(3,1);
errab = zeros(3,1);
errac = zeros(3,1);
errbb = zeros(3,1);
errbc = zeros(3,1);
errcc = zeros(3,1);

for i=1:1:row
    err(i) = data(i)-(a/(1+b*exp(-1*c+time(i))));
    erra = -1/(1+b*exp(-1*c*time(i)));
    errb = a*exp(-1*c*time(i))/((1+b*exp(-1*c*time(i)))^2);
    errc = -1*a*b*time(i)*exp(-1*c*time(i))/((1+b*exp(-1*c*time(i)))^2);
    errab = exp(-1*c*time(i))/((1+b*exp(c*time(1)))^2);
    errac = -b*time(i)*exp(-c*time(i))/((1+b*exp(-c*time(i)))^2);
    errbb = -2*a*exp(-c*time(i))^2/(1+b*(exp(-c*time(i))^3));
    errbc = (2*a*(exp(-c*time(i)))^2*b*time(i))/(1+b*exp(-c*time(i)))^3-(a*time(i)*exp(-c*time(i)))/(1+b*exp(-c*time(i)))^2;
    errcc = (-2*a*b^2*time(i)^2*(exp(-c*time(i)))^2)/(1+b*exp(-c*time(i)))^3+(a*b*time(i)^2*exp(-c*time(i)))/(1+b*exp(-c*time(i)))^2;
end

for i = 1:1:row
    df2(1,1) = df2(1,1)+2*(erra(i))^2;
    df2(1,2) = df2(1,2)+2*erra(i)*errb(i)+2*err(i)*errab(i);
    df2(1,3) = df2(1,3)+2*erra(i)*errc(i)+2*err(i)*errac(i);
    df2(2,2) = df2(2,2)+2*(errb(i))^2+2*err(i)*errbb(i);
    df2(2,3) = df2(2,3)+2*errb(i)*errc(i)+2*err(i)*errbc(i);
    df2(3,3) = df2(3,3)+2*(errc(i))^2+2*err(i)*errcc(i);
    % Hessian matrix is symetric because it doesn't matter which order you
    % take partial derivatives in
    df2(2,1) = df2(1,2);
    df2(3,2) = df2(2,3);
    df2(3,1) = df2(1,3);
end
end



