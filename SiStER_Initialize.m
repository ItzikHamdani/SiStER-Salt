% SiStER Initialize
% I. Hamdani (2017-2020) added:
%  	the option to run extra geometry features 
% 	add initial noise to plastic strain  
% 	track the sticky layer sediment interface

PARAMS.Nphase = Nphase; % for convenience

% construct staggered grids
[X,Y,x,y,xc,yc,dx,dy,Nx,Ny] = SiStER_initialize_grid(xsize,ysize,GRID);

% initialize marker arrays and positions
[xm, ym] = SiStER_initialize_marker_positions(xsize,ysize,dx,dy,Mquad);

% locate markers with respect to grid
[qd,icn,jcn] = SiStER_locate_markers_in_grid(xm,ym,x,y,dx,dy);

% initialize marker chain to track base of layer 1 (sticky layer)
Ntopo=PARAMS.Ntopo_markers;
topo_x=linspace(0,xsize,Ntopo);
topo_y=GEOM(1).bot*ones(size(topo_x));
topo_marker_spacing=mean(diff(topo_x)); % initial mean spacing of topography markers

% assign marker phases

[im] = SiStER_initialize_marker_phases(Nphase,GEOM,xm,ym);

% add extra geomtry features
if PARAMS.YNSaltGeometry
    Salt_extra_Geometry
end

xm_grd = X(1, 2:end-1);

% find seafloor level
b = repmat(topo_x, numel(xm_grd),1);
myArray = (b > xm_grd');
inter_i = sum( cumprod(myArray ==0, 2), 2) + 1;

dx0 = topo_x(inter_i) - xm_grd;
dx1 = xm_grd-topo_x(inter_i-1);
y_inter = ((topo_y(inter_i-1).*dx0) + (topo_y(inter_i).*dx1))...
            ./ (dx0+dx1);

[Y_INTER, ~] = meshgrid(y_inter, y);

% initialize marker plastic strain (to zero or noise) and strain rate (to one)
% c=0.01;
% ep=c*rand(size(xm)); NOISE

ep=zeros(size(xm));
epNH=ep;
epsIIm=ones(size(xm));

% initialize marker stresses
sxxm=zeros(size(xm));
sxym=sxxm;

% initialize marker index (a unique number to identify and track each marker)
idm=1:length(xm);

% initialize temperature structure on nodes
T=PARAMS.a0+PARAMS.a1*Y+PARAMS.a2*Y.^2+PARAMS.a3*Y.^3;
T(:, 2:end-1)=PARAMS.a0+PARAMS.a1*(Y(:, 2:end-1)-Y_INTER)+PARAMS.a2*Y(:, 2:end-1).^2+PARAMS.a3*Y(:, 2:end-1).^3;
T=T+PARAMS.amp*sin(2*pi*X/PARAMS.lam);
if PARAMS.ynTreset==1 % reset T=T0 in top layer
    T(T<PARAMS.T0)=PARAMS.T0;
end
% pass initial nodal T to markers
[Tm]=SiStER_interp_shear_nodes_to_markers(T,x,y,xm,ym,icn,jcn);
Tm0=Tm;

% initialize nodal strain rate and other useful arrays
eXX=zeros(size(X));
eXY=zeros(size(X));
EXX=zeros(size(X));
EXY=zeros(size(X));
vx=zeros(size(X));
vy=zeros(size(X));
v=vx;
p=1e12*ones(size(EXX));  %initialize to be high so plasticity doesnt activate at t=1, pit=1;
etan_new=zeros(Ny,Nx);
%-------------------------------------------------------------------------
% initialize dt_m small to keep things elastic & no plasticity at t=1, G.Ito
%-------------------------------------------------------------------------
if (exist('dt_m','var')==0)
    dt_m=1e2;
end


