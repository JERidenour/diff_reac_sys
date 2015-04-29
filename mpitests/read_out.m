%%
clear
%% Read info about matrix
minfo = dlmread('pPnN.txt');

p = minfo(1)
P = minfo(2)
n = minfo(3)
N = minfo(4)

%% Read M as a column of P blocks of nxn matrices indexed 0 1 2 ... P-1
M = dlmread('out.txt');

%% Assemble U as NxN matrix with pxp blocks indexed [ 0 1 2 3 ; 4 5 6 7; ...; P-4 P-3 P-2 P-1]

U = zeros(N,N);

for i = 1 : P
    m_in = (i*n)-(n-1);
    m_out = m_in + n-1;
    
    
    a_in_i = floor((i-1)/p)*n+1
    a_out_i = a_in_i+n-1
    a_in_j = mod(i+p-1,p)*n+1
    a_out_j = a_in_j+n-1
    
    U(a_in_i:a_out_i,a_in_j:a_out_j) = M(m_in:m_out,1:n);
end

%% plot

contourf(reshape(U,[N N]))