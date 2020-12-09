% plot ep
% run this script with the visualize loopto plot the plastic strain 

% plot
plotName = 'ep';
fastscatter(xm(im>1)/1e3,ym(im>1) /1e3,epNH(im>1),'markersize',2);
% fastscatter(xm(im>1)/1e3,ym(im>1) /1e3,ep(im>1),'markersize',2);
axis ij;
axis equal;
grid off;
c = colorbar;
ylabel(c, '${e}_{p}$', 'Interpreter','latex', 'FontSize', 16);
title([num2str(years(seconds(time))/1e3) 'ky']);
% set x, y, color limits
%     xlim([ 2.5 10]);
%     ylim([1 3]);
% caxis([0 0.005]);

if exist('toSave') == 1
% save as png if requested or show 
    if toSave
        fileName = [plotName '_' iStr '.png'];
        saveas(hFig, fileName);
    else
        set(hFig, 'visible', 'on');
        pause(PAUSE);
    end
end
