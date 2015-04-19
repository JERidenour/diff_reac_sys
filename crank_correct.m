%% params
clear


n = 81;
hx = 1/(n-1);    % i,j = [1...n]
hx2 = hx^2;

%ht = 1
%ht = hx
ht = 0.25;

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

U = ones(n,n);
V = zeros(n,n);

blopsize = 5;
ublop = zeros(blopsize)+0.5;
vblop = zeros(blopsize)+0.25;

in = floor(n/3)-floor(n/10);

U(in:in+size(ublop,1)-1,in:in+size(ublop,1)-1) = ublop;
V(in:in+size(vblop,1)-1,in:in+size(vblop,1)-1) = vblop;

in = floor(2*(n/3)+12);

U(in:in+size(ublop,1)-1,in:in+size(ublop,1)-1) = ublop;
V(in:in+size(vblop,1)-1,in:in+size(vblop,1)-1) = vblop;

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

u = reshape(U,[N 1]);
v = reshape(V,[N 1]);

u_new = u;
v_new = v;

%% I-T and I+T matrices (constant)
% 
ru = (Du/hx2)*(ht/2);
rv = (Dv/hx2)*(ht/2);

ImTu = sparse(eye(size(A)))-ru*A;
IpTu = sparse(eye(size(A)))+ru*A;

ImTv = sparse(eye(size(A)))-rv*A;
IpTv = sparse(eye(size(A)))+rv*A;


%% iterate

maxtime = 50;

for t = 1:maxtime
    
% RHS vectors (update every step)

RHSu = IpTu*u + ht*(-u.*(v.^2) + F*(1-u));
RHSv = IpTv*v + ht*(u.*(v.^2) - (F+k)*v);

u_new = ImTu\RHSu;
v_new = ImTv\RHSv;

u = u_new;
v = v_new;
    
% if mod(t,100)==0
%    contourf(reshape(u,[n n]))
%    pause(0.001)
% end



end
contourf(reshape(u,[n n]))