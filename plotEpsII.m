% plot strain rate
% A general script for plotting the second invariant of strain rate markers of a simulation iteration 
% I. Hamdani 2017

% plot
plotName = 'epsII';
% plot log values
fastscatter(xm(im>1)/1e3,ym(im>1) /1e3,log10(epsIIm(im>1)),'markersize',2);
% plot actual values
% fastscatter(xm(im>1)/1e3,ym(im>1) /1e3,(epsIIm(im>1)),'markersize',2);

axis ij;
axis equal;
grid off;
colorbar;
title([num2str(years(seconds(time))/1e3) 'ky']);

% set x, y, color limits
% xlim([40 50]);
% ylim([0 4]);
caxis([-18 -12]);

% save as png if requested or show on screen 
if exist('toSave') == 1
    if toSave
        fileName = [plotName '_' num2str(t) '.png'];
        saveas(gcf, fileName);
    else
        set(gcf, 'visible', 'on');
        pause(PAUSE);
    end
end