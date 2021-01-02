% plot ep
% A general script for plotting the plastic strain markers of a simulation iteration
% I. Hamdani 2017
 
% plot
plotName = 'ep';
% plot plastic strain w/o healing
fastscatter(xm(im>1)/1e3,ym(im>1) /1e3,epNH(im>1),'markersize',2);
% plot plastic strain w healing
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
