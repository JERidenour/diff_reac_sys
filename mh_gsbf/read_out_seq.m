clear

RU = dlmread('u_seq.txt');
RV = dlmread('v_seq.txt');
%%


figure(2)
subplot 121
contourf(RU)
axis equal
subplot 122
contourf(RV)
axis equal