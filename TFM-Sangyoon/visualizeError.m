function []=visualizeError(f,d,dispDetec,imgPath,method)

% make interpolated surface
% dispDetec_interp   = csapi({f,d},mean(dispDetec,3));
% figure, fnplt( dispDetec_interp )
% fu = f(1):f(end);
% du = d(1):.1:d(end);
% h1 = figure;
% dispDetec_interp_discrete = fnval(dispDetec_interp,{fu,du});
% pcolor(du,fu,dispDetec_interp_discrete)
% shading interp
% colormap hot
% hold on
% 
% contour(du,fu,dispDetec_interp_discrete,[1 1],'LineWidth',3,'LineColor',[32/255,77/255,2/255]);
% colorbar('location','eastoutside')

h1 = figure;
meanDispDetec = mean(dispDetec,3);
if strcmp(method,'pcolor')
    pcolor(d,f,meanDispDetec)
    shading interp
elseif strcmp(method,'pcolor_with_level1line')
    pcolor(d,f,meanDispDetec)
    shading interp
    hold on
    contour(d,f,meanDispDetec,[1 1],'LineWidth',3,'LineColor',[63/255,162/255,10/255]);
elseif strcmp(method, 'contourf')
    contourf(d,f,meanDispDetec,0:0.1:1);%[1 1],'LineWidth',3,'LineColor',[63/255,162/255,10/255]);    
else
    error('specify method with pcolor, pcolor_with_level1line or contourf')
end
dispDetecMax = max(meanDispDetec(:));
caxis([0 dispDetecMax])
colormap hot
colorbar('location','eastoutside')
xlabel('Adhesion size (pixel)')
ylabel('Input Force Magnitude (Pa)')

tifPath = [imgPath filesep 'tifs'];
figPath = [imgPath filesep 'figs'];
epsPath = [imgPath filesep 'eps'];
if ~exist(tifPath,'dir')
    mkdir(tifPath);
end
if ~exist(figPath,'dir')
    mkdir(figPath);
end
if ~exist(epsPath,'dir')
    mkdir(epsPath);
end

hgsave(h1,strcat(figPath,filesep,inputname(3)),'-v7.3')
hgexport(h1,strcat(tifPath,filesep,inputname(3)),hgexport('factorystyle'),'Format','tiff')
print(h1,strcat(epsPath,filesep,inputname(3),'.eps'),'-depsc2')
