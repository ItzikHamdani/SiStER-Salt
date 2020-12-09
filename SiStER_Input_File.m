% SiStER_Input_File

% PARAMETERS FOR SALT GEOMETRY
% angles of wedge
PARAMS.alpha = deg;
PARAMS.beta = deg;
% thidness of layers
PARAMS.salt_thickness = 1e3;
PARAMS.sed_thickness = thick;
% viscosities
PARAMS.SALT_VISC = 10^eta(1);
if eta(2) == 50
    PARAMS.SED_VISC = PARAMS.SALT_VISC*50;
else
    PARAMS.SED_VISC = 10^eta(2);
end
% PARAMS.WATER_VISC = 1e15; 
PARAMS.WATER_VISC =PARAMS.SALT_VISC*1e-2;
% sediment linear density profile- rho(d) = sed_rho(1)+d*sed_rho(2)
PARAMS.SED_RHO = [1877.6 0.195];
% special time points
PARAMS.TIMEPOINTS = [101 seconds(years([4.073e3 5e6]))]; % TODO
PARAMS.TPACTION = [1 2 0]; % TODO

PARAMS.wedge_offset = 5e3;
PARAMS.D = 0.01; % m
PARAMS.dy = 50;
PARAMS.y_nodes = 7;
% PARAMS.wedge_len = (PARAMS.salt_thickness)/ (tand(PARAMS.beta) - tand(PARAMS.alpha));
if PARAMS.alpha == PARAMS.beta
    PARAMS.wedge_len = 20e3;
else
    PARAMS.wedge_len = (PARAMS.salt_thickness - PARAMS.y_nodes*PARAMS.dy)/ (tand(PARAMS.beta) - tand(PARAMS.alpha));
end
PARAMS.wedge_h = tand(PARAMS.alpha)*PARAMS.wedge_len;
PARAMS.max_water_thickness = 1500;
% PARAMS.min_water_thickness = PARAMS.max_water_thickness - PARAMS.wedge_h;
PARAMS.min_water_thickness = 300;
% PARAMS.water_thickness = 100*ceil(PARAMS.wedge_h/100) + PARAMS.min_water_thickness; 
PARAMS.water_thickness = PARAMS.wedge_h + PARAMS.min_water_thickness; 
% PARAMS.water_thickness = PARAMS.max_water_thickness; 

% DURATION OF SIMULATION AND FREQUENCY OF OUTPUT %%%%%%%%%%%%%%%%%%%%%%%%%%
Nt=nt; % max number of time iterations
dt_out=1; % output files every "dt_out" iterations

% DOMAIN SIZE AND GRIDDING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xsize = PARAMS.wedge_offset + 3e3*ceil(PARAMS.wedge_len/1e3); % 2= 20 km of horizontal, 3 = 40 km etc
ysize= 100*(ceil((PARAMS.salt_thickness + PARAMS.sed_thickness + PARAMS.water_thickness)/100) +2);

% gridding- from 0 to GRID.x(1), grid size is GRID.dx(1)
% from GRID.x(1) to GRID.x(2), grid size is GRID.dx(1) etc...
% same for y
GRID.dx(1)=50;
GRID.x(1)=PARAMS.wedge_offset + 100*ceil(PARAMS.wedge_len/100) + PARAMS.wedge_offset;
GRID.dx(2)=100;
GRID.x(2)=xsize-10e3;
GRID.dx(3)=50;


GRID.dy(1)=50;
GRID.y(1)=900;
% GRID.dy(2)=25;
% GRID.y(2)=2.5e3;
% GRID.dy(3)=50;

% LAGRANGIAN MARKERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Mquad=4; % number of markers in the smallest quadrant
Mquad_crit=2; % minimum number of markers allowed in smallest quadrant (for reseeding)

% GEOMETRY %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Nphase=3; % number of phases

% phase 1 Water
GEOM(1).type=1; % 1 = layer (then specify top and bot) or 2 = circle % 1 = layer (then specify top and bot) or 2 = circle (then specify center and radius)
GEOM(1).top=0;
GEOM(1).bot=PARAMS.water_thickness;

% phase 2 Rock/Sediment
sed_phase_n=2;
GEOM(2).type=1; % 1 = layer (then specify top and bot) or 2 = circle % 1 = layer (then specify top and bot) or 2 = circle (then specify center and radius)
GEOM(2).top=GEOM(1).bot;
GEOM(2).bot=ysize;

% phase 3 Salt+wedge
wedge_phase_n=3;
GEOM(3).type=3; % 1 = layer (then specify top and bot) or 2 = circle (then specify center and radius)
GEOM(3).left = 5e3;
GEOM(3).right = xsize;
GEOM(3).top = PARAMS.sed_thickness + PARAMS.water_thickness;
GEOM(3).bot = GEOM(3).top + PARAMS.salt_thickness;
PARAMS.YNSaltGeometry = 1; % addition of wedges or other non rectangular shape

% MATERIAL PROPERTIES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% creep laws of the form: pre^(-1/n)*epsII^((1-n)/n)*exp(E/(nRT))
% harmonically averaging diffusion creep, dislocation creep 
% (and plastic creep to simulate brittle failure)

% phase 1 WATER/AIR
MAT(1).phase=1;
% density parameters
MAT(1).rho0=1000;   
MAT(1).alpha=0;
% thermal parameters
MAT(1).k=3;
MAT(1).cp=1000;
% elasticity 
MAT(1).G=1e18;
% diffusion creep parameters
MAT(1).pre_diff=.5/PARAMS.WATER_VISC;
MAT(1).Ediff=0;
MAT(1).ndiff=1;
% dislocation creep parameters
MAT(1).pre_disc=.5/PARAMS.WATER_VISC;
MAT(1).Edisc=0;
MAT(1).ndisc=1;
% plasticity
MAT(1).mu=0.6;
MAT(1).mumin=0.3;
MAT(1).Cmax=5e6;
MAT(1).Cmin=0.01e6;
MAT(1).ecrit=0.1;
MAT(1).YNPorePress=0;


% phase 2 ROCK
MAT(2).phase=2;
% density parameters
MAT(2).rho0=PARAMS.SED_RHO(1);
MAT(2).alpha=0;
% thermal parameters
MAT(2).k=3;
MAT(2).cp=1000;
% elasticity 
MAT(2).G=30e9;
% diffusion creep parameters
MAT(2).pre_diff=.5/PARAMS.SED_VISC ;
MAT(2).Ediff=0;
MAT(2).ndiff=1;
% dislocation creep parameters
MAT(2).pre_disc=.5/PARAMS.SED_VISC;
MAT(2).Edisc=0;
MAT(2).ndisc=1;
% plasticity
MAT(2).mu=0.53;
MAT(2).mumin=0.265;
MAT(2).Cmax=50e3;
MAT(2).Cmin=0.01e3;
MAT(2).ecrit=0.01;
MAT(2).YNPorePress=1;


% phase 3 SALT
MAT(3).phase=3;
% density parameters
MAT(3).rho0=2100;
MAT(3).alpha=0;
% thermal parameters
MAT(3).k=5;
MAT(3).cp=1000;
% elasticity 
MAT(3).G=30e9;
% diffusion creep parameters
MAT(3).pre_diff=.5/PARAMS.SALT_VISC;
MAT(3).Ediff=0;
MAT(3).ndiff=1;
% dislocation creep parameters
MAT(3).pre_disc=.5/PARAMS.SALT_VISC;
MAT(3).Edisc=0;
MAT(3).ndisc=1;
% plasticity
MAT(3).mu=0.6;
MAT(3).mumin=0.3;
MAT(3).Cmax=5e6; % TODO Cohesion 
MAT(3).Cmin=0.01e6;
MAT(3).ecrit=0.1;
MAT(3).YNPorePress=0;

% ADDITIONAL PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PARAMS.PATH = filepath;
PARAMS.YNThreshold = 0;
PARAMS.SLOPE_ANGLE = 0;
PARAMS.YNElast=0; % elasticity on (1) or off (0)
PARAMS.YNPlas=0; % plasticity on (1) or off (0)
PARAMS.YNSEDIMENT = 0; % sedimentation on (1) or off (0)
PARAMS.tau_heal=3e13; % healing time for plasticity (s) TODO
PARAMS.gx=9.8*sind(PARAMS.SLOPE_ANGLE); % gravity along x
PARAMS.gy=9.8*cosd(PARAMS.SLOPE_ANGLE); % gravity along y
PARAMS.fracCFL=0.05; % distance by which a marker is allowed to move over a time step, expressed as a fraction of the smallest cell size
PARAMS.R=8.314; % gas constant
PARAMS.etamax=PARAMS.SED_VISC*10; % maximum viscosity
% PARAMS.etamin=PARAMS.SALT_VISC*1e-3; % minimum viscosity PARAMS.WATER_VISC*0.1
PARAMS.etamin=PARAMS.WATER_VISC*1e-1; % minimum viscosity 
PARAMS.Tsolve=0; % yes (1) or no (0) solve for temperature
% initial temperature profile, polynomial with depth 
% T = a0 + a1*y+a2*y^2+a3*y^3+amp*sin(2*pi*X/lam)
% (make sure it matches the BCs)
PARAMS.a0=14;
PARAMS.a1=20/1000;
PARAMS.a2=0;
PARAMS.a3=0;
PARAMS.amp=0; % amplitude of sinusoidal perturbation
PARAMS.lam=1; % wavelength of sinusoidal perturbation
PARAMS.ynTreset=1; % if ==1, reset T=T0 where im==1 (sticky layer)
PARAMS.T0=50;
% reference values for the constant diffusivity thermal solver
% (kappa = kref / (rhoref*cpref))
PARAMS.rhoref=MAT(2).rho0; 
PARAMS.kref=3;
PARAMS.cpref=1000;

% TOPOGRAPHY EVOLUTION (interface between rock and sticky air/water layer)
PARAMS.Ntopo_markers=1000; % number of markers in marker chain tracking topography
PARAMS.YNSurfaceProcesses=0; % surface processes (diffusion of topography) on or off
PARAMS.topo_kappa=1e-8; % diffusivity of topography (m^2/s)


% Solver iterations
PARAMS.Npicard_min=20; % minimum number of Picard iterations per time step
PARAMS.Npicard_max=75; % maximum number of Picard iterations per time step
PARAMS.conv_crit_ResL2=1e-3;
PARAMS.pitswitch=0; % number of Picard iterations at which the solver switches to quasi-Newton


% BOUNDARY CONDITIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% pressure
PARAMS.p0cell=0; % pressure in the top-left corner of the domain (anchor point)


% flow

% boundary conditions
% entries in BC correspond to
% 1/ rollers? (1=yes, 0=no)
% 2/ type of velocity normal to boundary (0=constant)
% 3/ value of normal velocity 

BC.top=[1 0 0];
BC.bot=[0 0 0 0];
BC.left=[1 0 0];
BC.right=[1 0 0];
BC.bot_profile=zeros(1,xsize/min(GRID.dx)+1);
PARAMS.BalanceStickyLayer=1; % if set to 1, the code will reset the inflow 
% / outflow BCs to balance the inflow / outflow of sticky layer material,
% and rock separately, based on the position of the sticky layer / air
% interface TODO


% thermal 

% entries in BCtherm correspond to
% 1/ type? (1=Dirichlet, 2=Neumann)
% 2/ value
BCtherm.top=[1 0];
BCtherm.bot=[1 1000];
BCtherm.left=[2 0];
BCtherm.right=[2 0];
