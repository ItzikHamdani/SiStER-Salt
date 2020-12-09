function [rhom, lambda]=SiStER_get_density(im,Tm,MAT,ym,PARAMS, GEOM, xm, topo_x, topo_y)
% [rho]=SiStER_get_density(im,Tm,MAT)
% obtain density from temperature and material identity
% 
% edited by B. Klein, 9/20/2016 to remove struct indexing concatenating

T0=0;
rhom = zeros(size(im));
lambda = zeros(size(im));

types = unique(im);
for i = 1:length(types)
    logical = im == types(i);
    rho0 = MAT(types(i)).rho0;
    alpha = MAT(types(i)).alpha;
    rhom(logical) = rho0.*(1-alpha.*(Tm(logical)-T0));
end


% OLD WAY (slow)
%T0=0;
%rhom = [MAT(im).rho0].*(1-[MAT(im).alpha].*(Tm-T0));


% depth function - rho(y) = 1600 + y*10/18
wedge_phase_n = length(MAT);
im_rho = im~=1 & im~= wedge_phase_n;
xm_rho = xm(im_rho);

% find seafloor level
b = repmat(topo_x, numel(xm_rho),1);
myArray = (b > xm_rho');
inter_i = sum( cumprod(myArray ==0, 2), 2) + 1;

dx0 = topo_x(inter_i) - xm_rho;
dx1 = xm_rho-topo_x(inter_i-1);
y_inter = ((topo_y(inter_i-1).*dx0) + (topo_y(inter_i).*dx1))...
            ./ (dx0+dx1);

rhom(im_rho) = PARAMS.SED_RHO(1) + PARAMS.SED_RHO(2)*(ym(im_rho)-y_inter)/cosd(PARAMS.SLOPE_ANGLE);

if MAT(2).YNPorePress == 1
    % calculate only above salt layer
    im_pp = im_rho & ym < GEOM(wedge_phase_n).bot;
    % indices only for depth above salt bot
    is_pp = ym(im_pp) < GEOM(wedge_phase_n).bot;
    d1 = ym(im_pp)-y_inter(is_pp);

    % For Constant Density
    %           lambda(im_pp) = (ym(im_pp)*MAT(1).rho0) ./... 
    %                  (MAT(1).rho0*y_inter(is_pp) + d1(is_pp)*MAT(2).rho0);
    
    %     For Depth Dependent Density
    lambda(im_pp) = (ym(im_pp)*MAT(1).rho0) ./... 
             (MAT(1).rho0*y_inter(is_pp) + d1(is_pp).*(PARAMS.SED_RHO(1) + d1(is_pp)*PARAMS.SED_RHO(2)/2));        
end
	