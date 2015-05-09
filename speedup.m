%% sparse matrix version
T1 = 793.15;
T4 = 304.39;
T16 = 40.87;
T64 = 9.6;
T100 = 7.4;


P = [1 4 16 64 100];
Sp = [T1/T1 T1/T4 T1/T16 T1/T64 T1/T100];

%get a linear fit
coeffs = polyfit(P, Sp, 1);
% Get fitted values
fittedX = linspace(0, max(P));
fittedY = polyval(coeffs, fittedX);
% Plot the fitted line
hold on;
plot(fittedX, fittedY, 'r-');

plot(P,Sp, '*')
title('Parallel speedup for 400x400 gray-scott system, sparse matrix method')
xlabel('Number of processes')
ylabel('Speedup')
%axis([0 max(P)+10 0 max(P)+10])


%% kernel version
T1 = 113;
T4 = 94.26;
T16 = 9.778;
T64 = 4.248;
T100 = 5.556;


P = [1 4 16 64 100];
Sp = [T1/T1 T1/T4 T1/T16 T1/T64 T1/T100];

%get a linear fit
coeffs = polyfit(P, Sp, 2);
% Get fitted values
fittedX = linspace(0, max(P));
fittedY = polyval(coeffs, fittedX);
%alm = 1.284;
% Plot the fitted line
hold on;
plot(fittedX, fittedY, 'r-');
plot(P,Sp, '*')
%plot(fittedX, alm, 'k')
title('Parallel speedup for 400x400 gray-scott system, kernel method')
xlabel('Number of processes')
ylabel('Speedup')
%axis([0 max(P)+10 0 max(P)+10])