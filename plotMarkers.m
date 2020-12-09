% plot markers
% run this script with the visualize loopto plot the marker positions 

% plot
plotName = 'marker';
fastscatter(xm/1e3,ym/1e3,im,'markersize',2);
axis ij;
axis equal;
grid off;
% title(t);
colorbar;
title([years(seconds(time))/1e3 'ky']);
% set x, y, color limits
%     xlim([ 40 50]);
%     ylim([0 3]);
caxis([1 max(im)]);

% save as png if requested
if toSave
    fileName = [plotName '_' iStr '.png'];
    saveas(hFig, fileName);
else    
    set(hFig, 'visible', 'on');
    pause(PAUSE);
end
