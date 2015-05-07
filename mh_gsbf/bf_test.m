%% params
clear


n = 64;
hx = 1/(n-1)    % i,j = [1...n]
hx2 = hx^2;

ht = 0.25

Du = 2*1e-5;
Dv = 1e-5;
F = 0.026;
k = 0.0550;

%% sparse block tridiagonal matrix A for laplace with periodic boundaries

N = n^2;
% [main first first second second cornerBlock cornerBlock diagBlockCorner diagBlockCorner]
i = [1:N 2:N 1:N-1 n+1:N 1:N-n N-n+1:N 1:n n:n:N 1:n:N];    % rows
j = [1:N 1:N-1 2:N 1:N-n n+1:N 1:n N-n+1:N 1:n:N n:n:N];    % cols

subd1 = repmat([ones(1,n-1) 0],1,n);    % first
subd2 = repmat(ones(1,n),1,n-1);        % second
subd3 = ones(1,n);                      % corner block
%subd4 = [repmat([1 zeros(1,n-1)],1,n-1) 1]; % each diagonal block corner
subd4 = ones(1,n);                      % each diagonal block corner

s = [ (-4*ones(1,N)) subd1(1:end-1) subd1(1:end-1) subd2 subd2 subd3 subd3 subd4 subd4];  % vals
A = sparse(i,j,s);

%% sparse block tridiagonal matrix A for laplace no boundaries

% N = n^2;
% % [main first first second second cornerBlock cornerBlock diagBlockCorner diagBlockCorner]
% i = [1:N 2:N 1:N-1 n+1:N 1:N-n];    % rows
% j = [1:N 1:N-1 2:N 1:N-n n+1:N];    % cols
% 
% subd1 = repmat([ones(1,n-1) 0],1,n);    % first
% subd2 = repmat(ones(1,n),1,n-1);        % second
% %subd3 = ones(1,n);                      % corner block
% 
% s = [ (-4*ones(1,N)) subd1(1:end-1) subd1(1:end-1) subd2 subd2];  % vals
% A = sparse(i,j,s);

%% u = reshape(U) v = reshape(V)

% U = ones(n,n);
% V = zeros(n,n);
% 
% blopsize = 5;
% ublop = zeros(blopsize)+0.5;
% vblop = zeros(blopsize)+0.25;
% 
% in = floor(n/3)-floor(n/10);    % [15 19]
% 
% U(in:in+size(ublop,1)-1,in:in+size(ublop,1)-1) = ublop;
% V(in:in+size(vblop,1)-1,in:in+size(vblop,1)-1) = vblop;
% 
% in = floor(2*(n/3)+12); % [54 58]
% 
% U(in:in+size(ublop,1)-1,in:in+size(ublop,1)-1) = ublop;
% V(in:in+size(vblop,1)-1,in:in+size(vblop,1)-1) = vblop;
% 
% u = reshape(U,[N 1]);
% v = reshape(V,[N 1]);
% 
% u_new = u;
% v_new = v;

%%

U = ones(n,n);
V = zeros(n,n);

U(1:4,1:4) = 0.5;
V(1:4,1:4) = 0.25;
U(64-15:64-15+3,64-15:64-15+3) = 0.5;
V(64-15:64-15+3,64-15:64-15+3) = 0.25;

u = reshape(U,[N 1]);
v = reshape(V,[N 1]);

u_new = u;
v_new = v;
%%
% JR initial values
% U = ones(n);
% V = zeros(n);
% c = round(n/2-n/20);
% c2 = c+round(n/9);
% c3 = c-round(n/9);
% b = round(n/10);
% V(c3:c3+2*b, c2:c2+b) = 0.25;
% V(c2:c2+b, c3:c3+2*b) = 0.25;
% U(c3:c3+2*b, c2:c2+b) = 0.5;
% U(c2:c2+b, c3:c3+2*b) = 0.5;
%% add some small block

% U = ones(n,n);
% V = zeros(n,n);
% 
% U(4:5,4) = 0.5;
% V(4:5,4) = 0.25;
% U(5,5) = 0.5;
% V(5,5) = 0.25;
% U(12:13,12:13) = 0.5;
% V(12:13,12:13) = 0.25;
% 
% u = reshape(U,[N 1]);
% v = reshape(V,[N 1]);
% 
% u_new = u;
% v_new = v;



%% T

% Tu = (Du/hx2) * A;
% Tv = (Dv/hx2) * A;


%% Forward

%maxtime = 500000;
maxtime = 10000;

for t = 1:maxtime

% u_new = u + ht*(Tu*u -u.*(v.^2) + F*(1-u));
% v_new = v + ht*(Tv*v + u.*(v.^2) - (F+k)*v);

u_new = u + (ht*Du/hx2)*A*u + ht*(-u.*(v.^2) + F*(1-u));

v_new = v + (ht*Dv/hx2)*A*v + ht*(u.*(v.^2) - (F+k)*v);

u = u_new;
v = v_new;

if mod(t,30)==0
        subplot 121
        contourf(reshape(u,[n n]))
        axis equal
        subplot 122
        contourf(reshape(v,[n n]))
        axis equal
        pause(0.001)
    end

end

%%
figure(1)
subplot 121
contourf(reshape(u,[n n]))
axis equal
subplot 122
contourf(reshape(v,[n n]))
axis equal
pause(0.001)


