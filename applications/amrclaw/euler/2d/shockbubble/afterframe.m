colormap(jet)

axis([0 2.0 0 0.5])
axis equal; 
axis tight

fprintf('%10s : %12.4e\n','qmin',qmin);
fprintf('%10s : %12.4e\n','qmax',qmax);

caxis([0.1 2.81]);

showpatchborders
setpatchborderprops('linewidth',1);
% showgridlines(1:5);
% hidepatchborders;
set(gca,'fontsize',16);

tstr = sprintf('ForestClaw : t = %12.4f',t/2);
title(tstr,'fontsize',16);


if (PlotParallelPartitions==1)
    showpatchborders;
end

prt = false;
NoQuery = 0;
if (prt)
    MaxFrames = 31;
    fname = sprintf('sbamrmesh%2.2d.png',Frame);
    disp(fname);
    print('-dpng',fname);
end

shg;

clear afterframe
clear mapc2m
clear parallelpartitions