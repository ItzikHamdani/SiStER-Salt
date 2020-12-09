% plot Sxy
% run this script with the visualize loopto plot the shear stress 

% plot
plotName = 'sxy';
% fastscatter(xm(im>1)/1e3,ym(im>1) /1e3,(sxym(im>1)),'markersize',2);
fastscatter(xm/1e3,ym /1e3,(sxym),'markersize',2);
% fastscatter(xm(im>1 & xm <20e3)/1e3,ym(im>1 & xm <20e3) /1e3,(sxym(im>1 & xm <20e3)),'markersize',2);
axis ij;
axis equal;
grid off;
colorbar;
title([years(seconds(time))/1e3 'ky']);
% set x, y, color limits
%     xlim([ 2.5 10]);
%     ylim([1 3]);
% caxis([1 3]);

% save as png if requested
% if toSave
%     fileName = [plotName '_' iStr '.png'];
%     saveas(hFig, fileName);
% else
%     set(hFig, 'visible', 'on');
%     pause(PAUSE);
% end
