%% Numerical set up
clear all
%physics constants
Du = 2e-5;
Dv = 1e-5;
f = 0.0180; 
k = 0.0510;

%numerical parameters
N = 256*256; %total nr of elements
n = sqrt(N); %dimension of global mesh
h = 1/(n-1); %spacial step length
maxiter = 50000; %number of iterations
dt = 0.19; %temporal step length

%parallelization parameters
Psq = 9; %number of processes
P = sqrt(Psq); %dimension of process mesh
N_inner = N/Psq; %number of elements per process
n_inner = sqrt(N_inner); %dimension of local mesh

%sigma
su = Du/(h*h);
sv = Dv/(h*h);

%A
e = ones(n,1);
T=spdiags([e -4*e e], -1:1, n, n);
T(1,n) = 1;
T(n,1) = 1;
I = eye(n,n);
A = kron(I,T);
e = ones(n*n,2);
A = spdiags(e,[-n n],A);
A = spdiags(e, [-(n-1)*n (n-1)*n], A);

%preallocate
u = ones(N,1);
v = zeros(N,1);
u_new = u;
v_new = v;

% Stability parameters
% Cu = -2/max(eigs(su*A));
% display(Cu)
% Cv = -2/max(eigs(sv*A));
% display(Cv)
% display(dt)

%% run the simulation

%initial values
r = floor(n/20);
c = floor(n/2);
for i=c-r:c+r
    for j=c-r:c+r
        u(i+n*j) = 0.5;
        v(i+n*j) = 0.25;
    end
end

for i=0:maxiter

    unew = u + dt*su*A*u + dt*(-u.*v.^2 + f*(1-u));
    vnew = v + dt*sv*A*v + dt*( u.*v.^2 - (f+k)*v);
    
    u = unew;
    v = vnew;
    
    if mod(i,200)==0
        U=reshape(u,n,n);
        contourf(U);
        pause(0.05)
    end
    
end

U = reshape(u,n,n);
contourf(U)