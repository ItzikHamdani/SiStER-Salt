% plot markers
% A general script for plotting the markers of a simulation iteration
% I. Hamdani 2017

% plot
plotName = 'marker';
fastscatter(xm/1e3,ym/1e3,im,'markersize',2);

axis ij;
axis equal;
grid off;
colorbar;
title([years(seconds(time))/1e3 'ky']);

% set x, y, color limits
%     xlim([ 40 50]);
%     ylim([0 3]);
caxis([1 max(im)]);

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