% fdtd1d_step2.m
%
% Add a simple Gaussian soft source to Basic FDTD Engine

% INITIALIZE MATLAB
close all;
clc;
clear all;

% UNITS
meters      = 1;
centimeters = 1e-2 * meters;
millimeters = 1e-3 * meters;
inches      = 2.54 * centimeters;
feet        = 12 * inches;
seconds     = 1;
hertz       = 1/seconds;
kilohertz   = 1e3 * hertz;
megahertz   = 1e6 * hertz;
gigahertz   = 1e9 * hertz;

% CONSTANTS
c0 = 299792458 * meters/seconds;
e0 = 8.8541878176e-12 * 1/meters;
u0 = 1.2566370614e-6 * 1/meters;

% OPEN A FIGURE WINDOW
figure('Color','w');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% DASHBOARD
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% SOURCE PARAMETERS
fmax = 5.0 * gigahertz;  % max freq for simulation

% GRID PARAMETERS
nmax = 1;
NLAM = 10;  % grid resolution, use 20 to reduce dispersion
%NLAM = 20;  % grid resolution, use 20 to reduce dispersion
NBUFZ = [100 100];  % buffer before and after a device

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% COMPUTE OPTIMIZED GRID
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% NOMINAL RESOLUTION
lam0 = c0/fmax;
dz    = lam0/nmax/NLAM;

% COMPUTE GRID SIZE
Nz    = sum(NBUFZ) + 3;

% COMPUTE GRID AXIS
za = [0:Nz-1]*dz;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% BUILD DEVICE ON GRID
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% INITIALIZE MATERIALS TO FREE SPACE
ER = ones(1,Nz);
UR = ones(1,Nz);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% COMPUTE THE SOURCE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% COMPUTE TIME STEP (dt)
nbc = sqrt(UR(1)*ER(1));
dt = nbc*dz/(2*c0);

% COMPUTE SOURCE PARAMETERS
tau = 0.5/fmax;   % duration
t0 = 5*tau;       % offset

% COMPUTE NUMBER OF TIME STEPS
tprop = nmax*Nz*dz/c0;
t     = 2*t0 + 3*tprop;
STEPS = ceil(t/dt);

% COMPUTE THE SOURCE
t      = [0:STEPS-1]*dt;
nz_src = round(Nz/2);
Esrc = exp(-((t - t0)/tau).^2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% INITIALIZE FDTD PARAMETERS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% COMPUTE UPDATE COEFFICIENTS
mEy = (c0*dt)./ER;
mHx = (c0*dt)./UR;

% INITIALIZE FIELDS
Ey = zeros(1,Nz);
Hx = zeros(1,Nz);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PERFORM FDTD ANALYSIS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%
% MAIN LOOP -- ITERATIVE OVER TIME
%
for T = 1 : STEPS

    % UPDATE H FROM E
    for nz = 1 : Nz-1
       Hx(nz) = Hx(nz) + mHx(nz)*( Ey(nz+1) - Ey(nz) )/dz; 
    end
    Hx(Nz) = Hx(Nz) + mHx(Nz)*( 0 - Ey(Nz) )/dz; % Dirichlet BC
    
    % UPDATE E FROM H
    Ey(1) = Ey(1) + mEy(1)*( Hx(1) - 0 )/dz; % Dirichlet BC
    for nz = 2 : Nz
        Ey(nz) = Ey(nz) + mEy(nz)*( Hx(nz) - Hx(nz-1) )/dz;
    end
    
    % Inject E Source (Soft Source)
    Ey(nz_src) = Ey(nz_src) + Esrc(T);
    
    % Show Status
    if ~mod(T,1)
        
        % show fields
        draw1d(ER,Ey,Hx,dz);
        xlim([dz Nz*dz]);
        xlabel('z');
        title(['Field at step ' num2str(T) ' of ' num2str(STEPS)]);
        
        % draw graphics
        drawnow;
    end
end

    


