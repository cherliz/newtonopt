function [out] = logistic(x0)
[x,data]=runAll();
time = zeros(length(data),1);
for i=1:1:length(data)
    time(i) = i;
end
a=x0(1); b=x0(2); c=x0(3);
out = sum(data-(a/(1+a/b-1)*exp(-c*time)));
end

