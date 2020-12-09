% plot strain rate
% run this script with the visualize loopto plot the strain rate 

% plot


plotName = 'epsII';
fastscatter(xm(im>1)/1e3,ym(im>1) /1e3,log10(epsIIm(im>1)),'markersize',2);
% fastscatter(xm(im>1)/1e3,ym(im>1) /1e3,(epsIIm(im>1)),'markersize',2);
axis ij;
axis equal;
grid off;
colorbar;
title([num2str(years(seconds(time))/1e3) 'ky']);

% set x, y, color limits
% xlim([40 50]);
% ylim([0 4]);
%     xlim([ 2.5 10]);
%     ylim([1 3]);
caxis([-18 -12]);

% save as png if requested
if exist('toSave') == 1
    if toSave
        fileName = [plotName '_' iStr '.png'];
        saveas(hFig, fileName);
    else    
        set(hFig, 'visible', 'on');
        pause(PAUSE);
    end
end