%%
clear
%%
minfo = dlmread('pPnN.txt');

p = minfo(1)
P = minfo(2)
n = minfo(3)
N = minfo(4)

%%
M = dlmread('out.txt');


%% 

A = zeros(N,N);

% for i = 1 : p
%     
%    for j = 1 : p
%       in = (i*n)-(n-1);
%       out = in + n-1;
%       off = (j*n)-n
%       %[in out]
%       M(in+off:out+off,1:n)
%    end
%     
% end

for i = 1 : P
    m_in = (i*n)-(n-1);
    m_out = m_in + n-1;
    
    
    a_in_i = floor((i-1)/p)*n+1
    a_out_i = a_in_i+n-1
    a_in_j = mod(i+p-1,p)*n+1
    a_out_j = a_in_j+n-1
    
    A(a_in_i:a_out_i,a_in_j:a_out_j) = M(m_in:m_out,1:n);
end