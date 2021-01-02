% sister_timepoints
% this script is used to run commands/scripts based on specified timepoints if needed
% the cases codes should be given as an array in the SiStER_Input_File and
% and the corresponding "actions" should be called here.
% I. Hamdani 2017-2020

switch timepoint
    case 1 
        %PARAMS.SLOPE_ANGLE = 0.3;
        %PARAMS.gx=9.8*sind(PARAMS.SLOPE_ANGLE); % gravity along x
        %PARAMS.gy=9.8*cosd(PARAMS.SLOPE_ANGLE); % gravity along y
        
        %change sedRate  is = im ==1 & im < topo & im > line
%         ims = 2, zero params 
%         PARAMS.SED_RATE = 200/2.8e6/365/24/3600; % 200m/2.8Myr
        timeOffset = time;
        dt_out = 1;
        
        %PARAMS.YNElast=1; % elasticity on (1) or off (0)
        PARAMS.YNPlas=1; % plasticity on (1) or off (0)
    case 2
        PARAMS.YNPlas=1; % plasticity on (1) or off (0)
    case 3
        isSetTime = 0;
    case 4
        dt_out = 5;
end

