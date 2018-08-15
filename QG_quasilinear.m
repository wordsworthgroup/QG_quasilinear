%% QG_quasilinear.m
% reproduce the results of Wordsworth, Phys Fluids, 2009

close all
clear all

%% setup

% setup grid quantities
ny   = 2e2;  % number of spatial points []
nt   = 2e3; % number of temporal points []
Y    = 1.0; % size of domain in y direction [m]
T    = 5.0e3; % simulation time [s]
dt   = T/nt; % timestep [s]
dy   = Y/ny; % grid interval [m]
t_a  = (1:nt)*dt; % time vector [s]
y_a  = (1:ny)'*dy; % space vector [y]
iu   = 1:ny; % vector to access u from X
iQ   = ny+1:2*ny; % vector to access Q from X

% setup parameters
beta = 1.0; % gradient of planetary vorticity [1/ms]
k0   = 50;  % zonal wavenumber of eddies [1/m]
kapp = 0.0; % damping term [1/s]

% setup spatial derivative / integral matrices
e        = ones(ny,1);
I        = eye(ny);
dD       = spdiags([e -2*e e], -1:1, ny, ny);
dD(1,ny) = 1; 
dD(ny,1) = 1; % periodic boundary conditions
D2       = dD/dy^2; % double spatial derivative operator [1/m2]

% setup anonymous functions
gamm = @(u)    beta - D2*u; % gradient of zonal mean PV [1/ms]
Psi  = @(Q)    (D2 - k0^2*I)\Q; % eddy streamfunction variable [m2/s]
dQdt = @(u,Q) -1i*k0*(u.*Q + gamm(u).*Psi(Q)) - kapp*Q; % time rate of change of eddy PV variable [1/s/s]
dudt = @(u,Q) -(k0/2)*(imag(Psi(Q)).*real(Q) - imag(Q).*real(Psi(Q))) - kapp*u; % zonal acceleration [m/s/s]
dXdt = @(X)    [dudt(X(iu),X(iQ)); dQdt(X(iu),X(iQ))]; % time rate of change of state vector [m/s/s; 1/s/s]

% initial conditions
Q00    = 0.05; % initial eddy PV peak amplitude [1/s]
l0     = 6*pi; % initial eddy meridional wavenumber [1/m]
u0     = zeros(ny,1) + (rand(ny,1)-0.5)*1e-9; % initial zonal velocity [m/s]
Q0     = zeros(ny,1) + Q00*exp(1i*l0*y_a); % initial eddy PV variable [1/s]
u_hist = zeros(ny,nt); % zonal velocity history array [m/s]
Q_hist = zeros(ny,nt); % eddy PV variable history array [1/s]
X0     = [u0; Q0]; % initial state vector [m/s/s; 1/s/s]

%% get solution

fprintf('Calculating results...\n')
X = X0;
for n=1:nt 
    X           = RK4_fn(X,dXdt,dt);
    u_hist(:,n) = X(iu);
    Q_hist(:,n) = X(iQ);
end

%% display results

subplot(2,1,1)
surf(t_a,y_a,u_hist);
shading flat
view(2)
colorbar vert
xlabel('time [s]')
ylabel('y [m]')
title('u_{mean} [m/s]')

subplot(2,1,2)
surf(t_a,y_a,real(Q_hist));
shading flat
view(2)
colorbar vert
xlabel('time [s]')
ylabel('y [m]')
title('Re[Q] [1/s]')

colormap jet

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculate x_{n+1} using RK4

function xp = RK4_fn(x,f,h)

k1 = f(x);
k2 = f(x+h*k1/2);
k3 = f(x+h*k2/2);
k4 = f(x+h*k3);
xp = x + (h/6)*(k1 + 2*k2 + 2*k3 + k4);

end

