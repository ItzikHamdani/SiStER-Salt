% % Visualization of the run of a model specified by location

%% make videos
paths = {
         'D:\Salt model Results\20200223 simulation for stability analysis rev\1 deg\400m\18-22';
         'D:\Salt model Results\20200223 simulation for stability analysis rev\1 deg\600m\18-22';
         'D:\Salt model Results\20200223 simulation for stability analysis rev\1 deg\200m\16-19'
        };
n = length(paths);
% rate = ones(1, n)*4;
rate = 5;
% for i = 1:n
for i = 1
    cdir = dir(paths{i});
    filenames = {cdir(3:end).name};
    toPlot = filenames(contains(filenames, 'mat'));
    toPlot = str2double(regexprep(toPlot, '\D', ''));
    toPlot = sort(toPlot);
%     toPlot = 1:10;
%     video_epsII(paths{i}, toPlot, rate);
%     video_markers(paths{i}, toPlot, rate);
    video_markers_quiver(paths{i}, toPlot, rate);
%     video_sII(paths{i}, toPlot, rate);
%     video_ep(paths{i}, toPlot, rate);
%     video_sxx(paths{i}, toPlot, rate);
%     video_sxy(paths{i}, toPlot, rate);
end


%% 
% % optional - save figures as png files
% % use FILES to select what to plot, it is possible to plot several
% 
% myPath = 'D:\Salt model Results\20171224 low viscosity contrast';
% cd(myPath);

% 1- save plot png, 0 - show plot
toSave = 0;
PAUSE = 0.00;
% the iterations to plot
toPlot = 1:30:330;
paths = {'D:\Salt model Results\20181021\1500m water\with sed';
         'D:\Salt model Results\20181021\1500m water\with sed\sed xy';
         'D:\Salt model Results\20181021\2000m water\with sed';
         'D:\Salt model Results\20181021\2000m water\with sed\sed xy'};
% avg_v = zeros(1, numel(toPlot));
hFig = figure('visible', 'off', 'units', 'normalized', 'outerposition', [0 0 1 1]);
for i = 1:numel(toPlot)
    iStr = num2str(toPlot(i)); 
    load([paths{2} '\' iStr '.mat']);
    ky = num2str(time/365/24/3600/1000);
%     plotMarkers;
    plotEpsII;
%     plotSxy;
%     plotSxx;
%     plotEp;
%     avg_v(i) = mean(vx(20:58, 450));
%     figure(1)
%     contour(X,Y,p);
%     title(ky);
%     axis ij;
%     colorbar;
end
%% Velocity Field

clear; close all;
v = VideoWriter('velocity field.avi');
v.FrameRate = 7;
open(v);
set(gcf,'visible', 'off', 'units','normalized','position',[0 0 1 1]);
for j =[1 20:20:445]
    matFilename = sprintf('%d.mat', j);
    load(matFilename);
	quiver(X/1e3,Y/1e3,vx, vy, 'k','AutoScale','on', 'AutoScaleFactor', 15);
    axis ij;
    axis equal;
    grid off;
    xlim([0 15e3]);
    ylim([1e3 3e3]);
    title([num2str(time/365/24/3600/1000) 'kyrs']);
    
	frame=getframe(gcf); % leaving gcf out crops the frame in the movie.
	writeVideo(v,frame);
end
close(v);


%% pressure Field
% 
% clear; close all;
% v = VideoWriter('pressure_field.avi');
% v.FrameRate = 7;
% open(v);
% set(gcf,'visible', 'off', 'units','normalized','position',[0 0 1 1]);
% p_bot = [];
% pl_bot = [];
% p_salt = [];
% pl_salt = [];
% t_bot50 = [];
% % p_salt80 = [];
% % pl_salt80 = [];
% % t_bot80 =[];
% for j =[20:20:100]
%     matFilename = sprintf('%d.mat', j);
%     load(matFilename);
%     t_bot50 = [t_bot50 time];
%     p_salt = [p_salt p(58, 507)];
%     pl_salt = [pl_salt p(58, 2)];
% %     t_bot80 = [t_bot80 time];    
% %     p_salt80 = [p_salt80 p(58, 567)];
% %     pl_salt80 = [pl_salt80 p(58, 2)];
% % 	contour(X,Y,p);
% %     axis ij;
% %     colorbar;
% %     title([num2str(time/365/24/3600/1000) 'kyrs']);
% %     
% % 	frame=getframe(gcf); % leaving gcf out crops the frame in the movie.
% % 	writeVideo(v,frame);
% end
% % gp_salt80 = (p_salt80 - pl_salt80)/80e3;
% gp_salt50 = (p_salt - pl_salt)/50e3;
% 
% close(v);
% 
% %%
% fastscatter(xm(im>1)/1e3,ym(im>1) /1e3,log10(dsxxm(im>1)),'markersize',2);
% axis ij;
% axis equal;
% grid off;
% colorbar;
