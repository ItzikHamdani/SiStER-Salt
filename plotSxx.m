% plot sxx
% A general script for plotting the deviatoric stress markers of a simulation iteration
% I. Hamdani 2017

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