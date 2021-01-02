% This script is called by SiStER_initialize. and used to run special geometry for the salt
%	that is not supported by the original code
%	I. Hamdani 2017-2020

salt = GEOM(wedge_phase_n);
sed = GEOM(sed_phase_n);

% sediment below line base of salt

% salt wedge between top and base lines
if PARAMS.alpha == PARAMS.beta 
	im((im==3) &...
    ym > PARAMS.min_water_thickness + PARAMS.sed_thickness + PARAMS.salt_thickness + (xm-salt.left)*tand(PARAMS.beta)) = 2;

    im((im==2 | im==1) & xm > PARAMS.wedge_offset & xm < PARAMS.wedge_offset + PARAMS.wedge_len &...
    ym > PARAMS.min_water_thickness + PARAMS.sed_thickness + (xm-salt.left)*tand(PARAMS.alpha) & ...
    ym < PARAMS.min_water_thickness + PARAMS.salt_thickness + PARAMS.sed_thickness + (xm-salt.left)*tand(PARAMS.beta)) = 3;
else
    im((im==3) &...
    ym > PARAMS.min_water_thickness + PARAMS.sed_thickness + PARAMS.y_nodes*PARAMS.dy + (xm-salt.left)*tand(PARAMS.beta)) = 2;

    im((im==2 | im==1) & xm < PARAMS.wedge_offset + PARAMS.wedge_len & xm > PARAMS.wedge_offset &...
    ym > PARAMS.min_water_thickness + PARAMS.sed_thickness + (xm-salt.left)*tand(PARAMS.alpha) & ...
    ym < PARAMS.min_water_thickness + PARAMS.sed_thickness + PARAMS.y_nodes*PARAMS.dy + (xm-salt.left)*tand(PARAMS.beta)) = 3;
end

% sediment cover above wege and under minimum water thickness 
im( im==1 & ...
    ym > PARAMS.min_water_thickness & ...
    ym > PARAMS.min_water_thickness + (xm-salt.left)*tand(PARAMS.alpha)) = 2;

% update bathymetry to be minimum left of wedge
im_topo = topo_x <= salt.left;
topo_y(im_topo) = PARAMS.min_water_thickness;

% update bathymetry to follow line of sediments above wedge
im_topo = topo_x < PARAMS.wedge_offset + PARAMS.wedge_len & topo_x > salt.left;
topoX = topo_x(im_topo);
topo_y(im_topo) = PARAMS.min_water_thickness + (topoX-salt.left)*tand(PARAMS.alpha);
