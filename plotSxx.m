% plot sxx
% run this script with the visualize loopto plot the normal stress

% plot
plotName = 'sxx';
fastscatter(xm(im>1)/1e3,ym(im>1)/1e3,(sxxm(im>1)),'markersize',2);
axis ij;
axis equal;
grid off;
c=colorbar;
ylabel(c, '${\sigma}_{xx}$', 'Interpreter','latex', 'FontSize', 16);
title([num2str(years(seconds(time))/1e3) 'ky']);
% set x, y, color limits
%     xlim([ 2.5 10]);
%     ylim([1 3]);
caxis([-3e7 1e7]);

% save as png if requested or show 
if exist('toSave') == 1
    if toSave
        fileName = [plotName '_' iStr '.png'];
        saveas(hFig, fileName);
    else
        set(hFig, 'visible', 'on');
        pause(PAUSE);
    end
end