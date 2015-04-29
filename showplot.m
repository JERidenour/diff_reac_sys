clear all
load pro_0.txt -ascii
load pro_1.txt -ascii
load pro_2.txt -ascii
load pro_3.txt -ascii

A=[pro_2 pro_3; pro_0 pro_1];
contourf(A)

%load u_vals.txt -ascii

% N=8*8;
% n=sqrt(N);
% P = 2;
% Psq = P*P; %4
% inner_tot = N/Psq; %16
% inner_dim = sqrt(inner_tot); %4

%u = ones(N,1)*8;

% for i=1:inner_dim
%     
% %     %first row
% %     u(i) = u_vals(i);
% %     u(i+inner_dim) = u_vals(i + inner_tot);
% %     
% %     %second row
% %     u(i+2*inner_dim) = u_vals(i+inner_dim);
% %     u(i+3*inner_dim) = u_vals(i+inner_dim + inner_tot);
% %     
% %     %third row
% %     u(i+4*inner_dim) = u_vals(i + 2*inner_dim);
% %     u(i+5*inner_dim) = u_vals(i + 2*inner_dim + inner_tot);
% %     
% %     %fourth row
% %     u(i+6*inner_dim) = u_vals(i + 3*inner_dim);
% %     u(i+7*inner_dim) = u_vals(i + 3*inner_dim + inner_tot);
% %     
% %     %fifth row
% %     u(i+8*inner_dim) = u_vals(i+8*inner_dim);
% %     u(i+9*inner_dim) = u_vals(i+3*inner_tot);
% %     
% %     %sixth row
% %     u(i+10*inner_dim) = u_vals(i + inner_dim + 2*inner_tot);
% %     u(i+11*inner_dim) = u_vals(i + inner_dim + 3*inner_tot);
% %     
% %     %seventh row
% %     u(i+12*inner_dim) = u_vals(i + 2*inner_dim + 2*inner_tot);
% %     u(i+13*inner_dim) = u_vals(i + 2*inner_dim + 3*inner_tot);
% %     
% %     %eighth row
% %     u(i+14*inner_dim) = u_vals(i + 3*inner_dim + 2*inner_tot);
% %     u(i+15*inner_dim) = u_vals(i + 3*inner_dim + 3*inner_tot);
%     
% %     u(i) = u_vals(i);
% %     u(i + P) = u_vals(i + Psq);
% %     u(i + 2*P) = u_vals(i + 2*Psq);
% %     u(i + 3*P) = u_vals(i + P);
% %     u(i + 4*P) = u_vals(i + P + Psq);
% %     u(i + 5*P) = u_vals(i + P + 2*Psq);
% %     u(i + 6*P) = u_vals(i + 2*P);
% %     u(i + 7*P) = u_vals(i + 2*P + Psq);
% %     u(i + 8*P) = u_vals(i + 2*P + 2*Psq);
% %     
% %     u(i + 9*P) = u_vals(3*Psq + i);
% %     u(i + 10*P) = u_vals(3*Psq + i + Psq);
% %     u(i + 11*P) = u_vals(3*Psq + i + 2*Psq);
% %     u(i + 12*P) = u_vals(3*Psq + i + P);
% %     u(i + 13*P) = u_vals(3*Psq + i + P + Psq);
% %     u(i + 14*P) = u_vals(3*Psq + i + P + 2*Psq);
% %     u(i + 15*P) = u_vals(3*Psq + i + 2*P);
% %     u(i + 16*P) = u_vals(3*Psq + i + 2*P + Psq);
% %     u(i + 17*P) = u_vals(3*Psq + i + 2*P + 2*Psq);
%     
% %     u(i + 18*P) = u_vals(6*Psq + i);
% %     u(i + 19*P) = u_vals(6*Psq + i + Psq);
% %     u(i + 20*P) = u_vals(6*Psq + i + 2*Psq);
% %     u(i + 21*P) = u_vals(6*Psq + i + P);
% %     u(i + 22*P) = u_vals(6*Psq + i + P + Psq);
% %     u(i + 23*P) = u_vals(6*Psq + i + P + 2*Psq);
% %     u(i + 24*P) = u_vals(6*Psq + i + 2*P);
% %     u(i + 25*P) = u_vals(6*Psq + i + 2*P + Psq);
% %     u(i + 26*P) = u_vals(6*Psq + i + 2*P + 2*Psq);
%     
% end

%U = reshape(u,n,n);
%contourf(U);
