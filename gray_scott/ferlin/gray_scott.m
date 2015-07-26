%% Numerical set up
clear all
%physics constants
Du = 2e-6;
Dv = 1e-6;
f = 0.012; 
k = 0.050;

%numerical parameters
N = 400*400; %total nr of elements
n = sqrt(N); %dimension of global mesh
h = 1/(n-1); %spacial step length
maxiter = 40000; %number of iterations
dt = 0.75; %temporal step length

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
%lambdas = eigs(su*A);
%lambda_max = max(lambdas);
%Cu = -2/lambda_max;
%display(lambda_max)
%display(Cu)
% Cv = -2/max(eigs(sv*A));
% display(Cv)
%display(dt)

%% run the simulation
str = sprintf('gray_scott_F%d_K%d.avi', f, k);
writerObj = VideoWriter(str);
writerObj.FrameRate = 50;
open(writerObj);
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
    
    if mod(i,40)==0
        U=reshape(u,n,n);
        contourf(U);
        frame = getframe;
        writeVideo(writerObj, frame);
        %pause(0.05)
   end
    
end
close(writerObj);
% x = 0:h:1;
% y = 0:h:1;
% U = reshape(u,n,n);
% contourf(x,y,U)
% str = sprintf('Gray-Scott system, %d iterations', maxiter);
% title(str)