function printPNGEPSFIG(fhandle,outputDirPlot,filename)
set(0,'CurrentFigure',fhandle);
print([outputDirPlot filesep filename '.png'],'-dpng');
print([outputDirPlot filesep filename '.eps'],'-depsc');
saveas(gcf,[outputDirPlot filesep filename '.fig'])
