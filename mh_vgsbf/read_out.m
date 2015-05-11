%%
clear
%% Read info about matrix
minfo = dlmread('pPnN.txt');

p = minfo(1)
P = minfo(2)
n = minfo(3)
N = minfo(4)

%% Read M as a column of P blocks of nxn matrices indexed 0 1 2 ... P-1
RU = dlmread('u.txt');
RV = dlmread('v.txt');

%% Assemble U as NxN matrix with pxp blocks indexed [ 0 1 2 3 ; 4 5 6 7; ...; P-4 P-3 P-2 P-1]

UMPI = zeros(N,N);
VMPI = zeros(N,N);

for i = 1 : P
    m_in = (i*n)-(n-1);
    m_out = m_in + n-1;
    
    
    a_in_i = floor((i-1)/p)*n+1
    a_out_i = a_in_i+n-1
    a_in_j = mod(i+p-1,p)*n+1
    a_out_j = a_in_j+n-1
    
    UMPI(a_in_i:a_out_i,a_in_j:a_out_j) = RU(m_in:m_out,1:n);
    VMPI(a_in_i:a_out_i,a_in_j:a_out_j) = RV(m_in:m_out,1:n);
end

%% plot
figure(2)
subplot 121
% contourf(reshape(UMPI,[N N]))

contourf(circshift(circshift(reshape(UMPI,[N N]),N/2,1 ),N/2,2))
title('U')
axis equal
subplot 122
% contourf(reshape(VMPI,[N N]))
contourf(circshift(circshift(reshape(VMPI,[N N]),N/2,1 ),N/2,2))
title('V')
axis equal

%%
