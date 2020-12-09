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
y_sed = 0; % TODO vectorize % move to if statements
x_sed = 0;
max_y_sed = min(ym(im==2));

choosePoints 

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
isSetTime = 1;
dt_m1 = 500*3600*24*365;

R_MAT = nan(Nt, PARAMS.Npicard_min);  % itzik TODO ASSES CONVERGENCE

while t <= Nt && time < PARAMS.TIMEPOINTS(end)% time loop
    disp(['STARTING ITERATION: ' num2str(t) ' out of ' num2str(Nt)])
    if PARAMS.YNSEDIMENT == 1
        y_sed = y_sed + dt_m*PARAMS.SED_RATE;
        x_sed = x_sed + dt_m*PARAMS.SED_RATEx;
        % TODO save indices of sedimentation and re-initialize marks strairate
        % and stresses.
        agg = ym > GEOM(wedge_phase_n).top - y_sed;
        prog = ym > max_y_sed & ym > tand(WEDGE1_SLOPE)*(xm-salt.left-x_sed) + salt.top;
        to_sediment = im ==1 & (agg | prog);
        im(to_sediment) = 2;
        epsIIm(to_sediment) = 0;
        sxxm(to_sediment) = 0;
        sxym(to_sediment) = 0;
        dt_m = 1e2; % TODO is necessary?
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
    
    timepoint = PARAMS.TPACTION(time-dt_m<PARAMS.TIMEPOINTS & time>PARAMS.TIMEPOINTS);
    % OUTPUT VARIABLES OF INTEREST (prior to rotation & advection)
    if (mod(t,dt_out)==0 && dt_out>0) || t==1 || t==Nt || ~isempty(timepoint) % SAVING SELECTED OUTPUT
        disp('SAVING SELECTED VARIABLES TO OUTPUT FILE') 
        filename=[PARAMS.PATH '\' num2str(t)];
        [etam]=SiStER_interp_shear_nodes_to_markers(etas,x,y,xm,ym,icn,jcn); % to visualize viscosity on markers
        %
        
        
        % save all vars for continuing simulation. optional - add condition
        % on time if want to stop at specific time
        % could save the find result, and keep additional array with
        % commands what to do in each of the timepoints - save, add tilt,
        % add sediments or whatever with a switch case code 
        % TODO timepoint = PARAMS.TPACTION(find(....))
%         timepoint = find(time-dt_m<PARAMS.TIMEPOINTS & time>PARAMS.TIMEPOINTS, 1);
%         if ~isempty(timepoint) 
%             save(filename);
%             sister_timepoints;
%         elseif t==Nt || t==1 
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

    