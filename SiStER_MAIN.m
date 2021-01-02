% SiStER_MAIN.m
%
% Simple Stokes solver with Exotic Rheologies
%
% Main routine doing initialization, time loop and outputs
%
%
% J.-A. Olive, B.Z. Klein, E. Mittelstaedt, M. Behn, G. Ito, S. Howell
% jaolive <at> ldeo.columbia.edu
% March 2011 - April 2017
% I. Hamdani (2017-2020) addded an option to:
% 		run specified scripts add specified timepoints, 
%		track selected variables also between saved model iterations outputs
%		add sedimentation
%		continue a simulation from an iteration (where the entire workspace was saved)


close all

% INITIALIZATION

% Input File: loads parameter values, model geometry, boundary conditions
if exist('running_from_SiStER_RUN','var')==0
    clear 
    InpFil = input('Input file ? ','s');
end
cur_iter = 1;
if exist('running_from_continue','var') && running_from_continue == 1
    load(InpFil);
    cur_iter = t+1;
    Nt = Nt + iters;
	PARAMS.TIMEPOINTS = [PARAMS.TIMEPOINTS, PARAMS.TIMEPOINTS(end)*2];

elseif running_from_SiStER_RUN == 1    
    [filepath,name,ext] = fileparts(InpFil);
    run(InpFil)

        % construct grid and initialize marker / node arrays
    SiStER_Initialize

    % BEGIN TIME LOOP %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    time=0;
end
t = cur_iter;

% initialize sedimentation thickness if sedimentation is on
if PARAMS.YNSEDIMENT == 1
	y_sed = 0; 
	x_sed = 0;
	max_y_sed = min(ym(im==2));
end

%%% use as an additional tracking of specific parameters between saved iterations
choosePoints % script that defines nPoints in the 2D grid to track, xi, yi_salt and y_sed are arrays of indices in the grid

salt_i = sub2ind(size(X),yi_salt, xi);
sed_i = sub2ind(size(X),yi_sed, xi);

sxx_ref_sed = nan(Nt, nPoints);
sxy_ref_sed = nan(Nt, nPoints);
EXX_ref_sed = nan(Nt, nPoints);
EXY_ref_sed = nan(Nt, nPoints);
vx_ref_sed = nan(Nt, nPoints);
vy_ref_sed = nan(Nt, nPoints);
t_ref = nan(Nt, 1);

sxx_ref_salt = nan(Nt, nPoints);
sxy_ref_salt = nan(Nt, nPoints);
EXX_ref_salt = nan(Nt, nPoints);
EXY_ref_salt = nan(Nt, nPoints);
vx_ref_salt = nan(Nt, nPoints);
vy_ref_salt = nan(Nt, nPoints);

% track CONVERGENCE
R_MAT = nan(Nt, PARAMS.Npicard_min);  

isSetTime = 1;
dt_m1 = 500*3600*24*365;

while t <= Nt && time < PARAMS.TIMEPOINTS(end)% time loop
    disp(['STARTING ITERATION: ' num2str(t) ' out of ' num2str(Nt)])
    
	% turn water/sticky layer to sediment - if sedimentation is on
	if PARAMS.YNSEDIMENT == 1
        % update the sticky layer-sediment interface location
		y_sed = y_sed + dt_m*PARAMS.SED_RATE;
        x_sed = x_sed + dt_m*PARAMS.SED_RATEx;
        % save indices of sedimentation and re-initialize markers strain rate
        % and stresses to appropriate values
        agg = ym > GEOM(wedge_phase_n).top - y_sed;
        prog = ym > max_y_sed & ym > tand(WEDGE1_SLOPE)*(xm-salt.left-x_sed) + salt.top;
        to_sediment = im ==1 & (agg | prog);
        im(to_sediment) = 2;
        epsIIm(to_sediment) = 0; % 
        sxxm(to_sediment) = 0;
        sxym(to_sediment) = 0;
    end    
    % update time
    time=time+dt_m;
    
    % Here we prepare nodal arrays to feed the Stokes solver 
    SiStER_material_props_on_nodes

    %%% SOLVE STOKES WITH NON-LINEAR RHEOLOGY HERE 
    SiStER_flow_solve
    
    % GET STRAIN RATE FROM CURRENT SOLUTION
    epsIIm=SiStER_interp_shear_nodes_to_markers(epsII_s,x,y,xm,ym,icn,jcn);
    
	[n2interp] = SiStER_interp_markers_to_shear_nodes(xm,ym,icn,jcn,qd,x,y,sxym);
    sxy = n2interp(1).data;  

    [n2interp] = SiStER_interp_markers_to_normal_nodes(xm,ym,icn,jcn,x,y,sxxm);
    sxx = n2interp(1).data; 
    
	% reference values:
    t_ref(t) = time;
    sxx_ref_sed(t, :) = sxx(sed_i);
    sxy_ref_sed(t, :) = sxy(sed_i);
    EXX_ref_sed(t, :) = EXX(sed_i);
    EXY_ref_sed(t, :) = EXY(sed_i);
    vx_ref_sed(t, :) = vx(sed_i);
    vy_ref_sed(t, :) = vy(sed_i);

    sxx_ref_salt(t, :) = sxx(salt_i);
    sxy_ref_salt(t, :) = sxy(salt_i);
    EXX_ref_salt(t, :) = EXX(salt_i);
    EXY_ref_salt(t, :) = EXY(salt_i);
    vx_ref_salt(t, :) = vx(salt_i);
    vy_ref_salt(t, :) = vy(salt_i);
    
	% find if this is a timepoint designated for running some script (e.g. change some values\parameters)
    timepoint = PARAMS.TPACTION(time-dt_m<PARAMS.TIMEPOINTS & time>PARAMS.TIMEPOINTS);
    % OUTPUT VARIABLES OF INTEREST (prior to rotation & advection)
    if (mod(t,dt_out)==0 && dt_out>0) || t==1 || t==Nt || ~isempty(timepoint) % SAVING SELECTED OUTPUT
        disp('SAVING SELECTED VARIABLES TO OUTPUT FILE') 
        filename=[PARAMS.PATH '\' num2str(t)];
        [etam]=SiStER_interp_shear_nodes_to_markers(etas,x,y,xm,ym,icn,jcn); % to visualize viscosity on markers
        %
        
        
        % save all vars for continuing simulation. 
        if t==Nt || t==1 
            save(filename);
        else
            save(filename,'X','Y','vx','vy','p','time','xm','ym','etam',...
                'rhom','BC','etan','Tm','im','idm','epsIIm','sxxm','sxym'...
                ,'ep','epNH','icn','jcn','qd','topo_x','topo_y','dt_m',...
                'Zn', 'Zs', 'sxxOLD', 'sxxOLD_s', 'sxyOLD', 'sxyOLD_n',...
                'dsxxm', 'dsxym', 'sxx', 'sxy', 'EXY', 'EXYOLD', 'EXX',...
                'EXXOLD', 'etas', 'eXX', 'eXY', 'R_MAT', ...
                'sxx_ref_sed', 'sxy_ref_sed', 'EXX_ref_sed', 'EXY_ref_sed', ...
                'vx_ref_sed', 'vy_ref_sed', 'sxx_ref_salt', 'sxy_ref_salt', ...
                'EXX_ref_salt', 'EXY_ref_salt', 'vx_ref_salt', ...
                'vy_ref_salt', 't_ref');
        end
    end
	
    % USE STRAIN RATE TO UPDATE STRESSES ON MARKERS
    SiStER_update_marker_stresses;
    
    % BUILD UP PLASTIC STRAIN IN YIELDING AREAS IF PLASTICITY IS ACTIVATED
    if (PARAMS.YNPlas==1) 
        SiStER_update_ep;
    end
    
    % SET ADVECTION TIME STEP BASED ON CURRENT FLOW SOLUTION
    
    [dt_m]=SiStER_set_timestep(dx,dy,vx,vy,PARAMS);
    
    if isSetTime
        dt_m = min([dt_m1, dt_m]);
    end
    

    % ROTATE ELASTIC STRESSES IN CURRENT FLOW FIELD
    if (PARAMS.YNElast==1) 
        SiStER_rotate_stresses;
    end
    
    % EVOLVE TEMPERATURE FIELD THROUGH DIFFUSION
    if PARAMS.Tsolve==1
        SiStER_thermal_update;
    end

    % MARKER ADVECTION, REMOVAL, AND ADDITION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    SiStER_move_remove_and_reseed_markers;
    % advect markers in current flow field
    % remove markers if necessary
    % add markers if necessary
    SiStER_update_topography_markers;
    % here we do the same for the marker chain that keeps track of topography
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
	
	%  run a script for this timepoint
	if ~isempty(timepoint) 
            save(filename);
            sister_timepoints;
    end
    disp('---------------')
    disp(['END OF ITERATION: ' num2str(t) ' out of ' num2str(Nt) ' - SIMULATION TIME: ' num2str(time/365.25/24/3600/1000) ' kyrs.'])
    disp('--------------------------------')
    disp('--------------------------------')
    
    t = t+1;
end

disp('FIN')

    