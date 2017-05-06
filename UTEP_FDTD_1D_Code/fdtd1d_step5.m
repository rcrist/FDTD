% fdtd1d_step5.m
%
% Add reflectance and transmittance calculation

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
NFREQ = 1000;
FREQ = linspace(0,fmax,NFREQ);

% GRID PARAMETERS
nmax = 1;
NLAM = 20;  % grid resolution, use 20 to reduce dispersion
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
s      = dz/(2*c0) + dt/2;
nz_src = 2;
Esrc   = exp(-((t - t0)/tau).^2);
A      = - sqrt(ER(nz_src)/UR(nz_src));
Hsrc   = A*exp(-((t - t0 + s)/tau).^2);

% plot(t,Esrc,'b');
% hold on;
% plot(t,-Hsrc,'-r')
% hold off;
% return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% INITIALIZE FDTD PARAMETERS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% INITIALIZE FOURIER TRANSFORMS
K   = exp(-1i*2*pi*dt*FREQ);  % Fourier transform kernal
REF = zeros(1,NFREQ);
TRN = zeros(1,NFREQ);
SRC = zeros(1,NFREQ);

% COMPUTE UPDATE COEFFICIENTS
mEy = (c0*dt)./ER;
mHx = (c0*dt)./UR;

% INITIALIZE FIELDS
Ey = zeros(1,Nz);
Hx = zeros(1,Nz);

% INITIALIZE BOUNDARY TERMS
H1=0; H2=0; H3=0;
E1=0; E2=0; E3=0;

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
    Hx(Nz) = Hx(Nz) + mHx(Nz)*( E3 - Ey(Nz) )/dz; % PABC
    
    % H Source
    Hx(nz_src-1) = Hx(nz_src-1) - mHx(nz_src-1)*Esrc(T)/dz;  %TF/SF correction
    
    % Record H-Field at Boundary
    H3 = H2;
    H2 = H1;
    H1 = Hx(1);
    
    % UPDATE E FROM H
    Ey(1) = Ey(1) + mEy(1)*( Hx(1) - H3 )/dz; % PABC
    for nz = 2 : Nz
        Ey(nz) = Ey(nz) + mEy(nz)*( Hx(nz) - Hx(nz-1) )/dz;
    end
    
    % E Source
    Ey(nz_src) = Ey(nz_src) - mEy(nz_src)*Hsrc(T)/dz;
    
    % Record E-Field at Bounday
    E3 = E2;
    E2 = E1;
    E1 = Ey(Nz);
    
    % Update Fourier Transforms
    for nf = 1 : NFREQ
       REF(nf) = REF(nf) + (K(nf)^T)*Ey(1);  % Intergrate over time
       TRN(nf) = TRN(nf) + (K(nf)^T)*Ey(Nz);
       SRC(nf) = SRC(nf) + (K(nf)^T)*Esrc(T);
    end
        
    % Show Status
    if ~mod(T,100)
        
        % show fields
        subplot(211);
        draw1d(ER,Ey,Hx,dz);
        xlim([dz Nz*dz]);
        xlabel('z');
        title(['Field at step ' num2str(T) ' of ' num2str(STEPS)]);
        
        % temporary normilization 
        R = abs(REF./SRC).^2;
        T = abs(TRN./SRC).^2;
        
        % show REF and TRN
        subplot(212);
        plot(FREQ,R,'-r');
        hold on;
        plot(FREQ,T,'b');
        plot(FREQ,R+T,':k');
        xlim([FREQ(1) FREQ(NFREQ)]);
        ylim([-0.1 1.5]);
        xlabel('Frequency (Hz)');
        title('REFLECTANCE AND TRANSMITTANCE');
        hold off;
        % draw graphics
        drawnow;
    end
end

% NORMALIZE FOURIER TRANSFORMS TO SOURCE
REF = abs(REF./SRC).^2;
TRN = abs(TRN./SRC).^2;
CON = REF + TRN;  % Check conservation of power

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% GENERATE A PROFESSIONAL LOOKING PLOT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% INITIALIZE FIGURE
clf;

% PLOT LINES
h = plot(FREQ/gigahertz,100*REF,'-r','LineWidth',2);
hold on;
plot(FREQ/gigahertz,100*TRN,'-b','LineWidth',2);
plot(FREQ/gigahertz,100*CON,':k','LineWidth',2);
hold off;

% SET VIEW
xlim([FREQ(1) FREQ(NFREQ)]/gigahertz);
ylim([0 105]);
h2 = get(h,'Parent');
set(h2,'FontSize',14,'LineWidth',2);
h = legend('Reflectance','Transmittance','Conservation');
set(h,'Location','NorthOutside');

% LABEL PLOT
xlabel('Frequency (GHz)');
ylabel('%', 'Rotation', 0, 'HorizontalAlignment','right');




