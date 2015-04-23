clear all
load u_vals.txt -ascii

N=81;
n=sqrt(N);
Psq = 9;
P = 3;
u = zeros(N,1);

for i=1:P
    
    u(i) = u_vals(i);
    u(i + P) = u_vals(i + Psq);
    u(i + 2*P) = u_vals(i + 2*Psq);
    u(i + 3*P) = u_vals(i + P);
    u(i + 4*P) = u_vals(i + P + Psq);
    u(i + 5*P) = u_vals(i + P + 2*Psq);
    u(i + 6*P) = u_vals(i + 2*P);
    u(i + 7*P) = u_vals(i + 2*P + Psq);
    u(i + 8*P) = u_vals(i + 2*P + 2*Psq);
    
    u(i + 9*P) = u_vals(3*Psq + i);
    u(i + 10*P) = u_vals(3*Psq + i + Psq);
    u(i + 11*P) = u_vals(3*Psq + i + 2*Psq);
    u(i + 12*P) = u_vals(3*Psq + i + P);
    u(i + 13*P) = u_vals(3*Psq + i + P + Psq);
    u(i + 14*P) = u_vals(3*Psq + i + P + 2*Psq);
    u(i + 15*P) = u_vals(3*Psq + i + 2*P);
    u(i + 16*P) = u_vals(3*Psq + i + 2*P + Psq);
    u(i + 17*P) = u_vals(3*Psq + i + 2*P + 2*Psq);
    
    u(i + 18*P) = u_vals(6*Psq + i);
    u(i + 19*P) = u_vals(6*Psq + i + Psq);
    u(i + 20*P) = u_vals(6*Psq + i + 2*Psq);
    u(i + 21*P) = u_vals(6*Psq + i + P);
    u(i + 22*P) = u_vals(6*Psq + i + P + Psq);
    u(i + 23*P) = u_vals(6*Psq + i + P + 2*Psq);
    u(i + 24*P) = u_vals(6*Psq + i + 2*P);
    u(i + 25*P) = u_vals(6*Psq + i + 2*P + Psq);
    u(i + 26*P) = u_vals(6*Psq + i + 2*P + 2*Psq);
    
end
U = reshape(u,n,n);
contourf(U);
