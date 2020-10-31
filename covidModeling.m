function [] = covidModeling()
x0=[100,1000000,100000000];
[x,data]=runAll();
time = zeros(length(data),1);
for i=1:1:length(data)
    time(i) = i;
end

[xbest,fbest,itrcount,status] = NewtonOpt(@logistic,@logisticDF,@logisticDF2,[400000;1000;.035],1,1,100,2);
xbest
%[9000000;105;0.025]
% Assemble the Logistic model based on the coefficients from NewtonOpt
f = @(x) xbest(1) ./ (1 + xbest(2).*exp(-xbest(3).*x));
trend = time(1)+1:0.1:time(end)+length(time);
% 60 time Prediction
mo2 = time(end):0.1:time(end)+60;

% Big thank you to Cory Snyder for showing me how to make cool plots!!!
figure(2);
hold on;
grid on;
plot(time,data, '-o');
plot(trend,f(trend),'LineWidth',1.5);
plot(mo2,f(mo2),'--','LineWidth',1.5);
title('Logistic Fit: Total Confirmed Cases in the US');
xlabel('Days');
ylabel('Confirmed Cases');
legend('Original Data','Logistic Curve','2 Month Forecast');
hold off;

end

