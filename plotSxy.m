% plot Sxy
% A general script for plotting the shear stress markers of a simulation iteration
% I. Hamdani 2017

% plot
plotName = 'sxy';
% plot w/o sticky layer
% fastscatter(xm(im>1)/1e3,ym(im>1) /1e3,(sxym(im>1)),'markersize',2);
% plot w sticky layer
fastscatter(xm/1e3,ym /1e3,(sxym),'markersize',2);

axis ij;
axis equal;
grid off;
colorbar;
title([years(seconds(time))/1e3 'ky']);

% set x, y, color limits
%     xlim([ 2.5 10]);
%     ylim([1 3]);
% caxis([1 3]);

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