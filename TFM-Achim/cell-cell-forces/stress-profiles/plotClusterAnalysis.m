function plotClusterAnalysis(constrForceField,forceField,imageFileList,frame,plotAll)

if nargin < 1 || isempty(constrForceField)
   [filename, pathname] = uigetfile({'*.mat';'*.*'}, ...
       'Select cellCellForces.mat to be used');
       fileStruct=load([pathname filesep filename]);
       constrForceField=fileStruct.constrForceField;
end

if nargin < 2 || isempty(forceField)
   [filename, pathname] = uigetfile({'*.mat';'*.*'}, ...
       'Select forceField.mat to be used');
       fileStruct=load([pathname filesep filename]);
       forceField=fileStruct.forceField;
end

%read in stack of image files:
if nargin < 3 || isempty(imageFileList)
   [filename, pathname] = uigetfile({'*.tif';'*.TIF';'*.*'}, ...
       'Select the first image file to be overlayed');   
   if ~ischar(filename) || ~ischar(pathname)
       return;
   end   
   imageFileList = getFileStackNames([pathname filesep filename]);
else
    isValid = 1;
    for i = 1:numel(imageFileList)
        isValid = isValid && exist(imageFileList{i}, 'file');
    end
    if ~isValid
        error('Invalid input files.');
    end
end

if nargin < 4 || isempty(frame)
    frame=input('Which frame should be analyzed: ');
    if isempty(frame)
        return
    end
end

% first test if the cluster analysis has been performed on this frame. For
% example, there is only a single cell:
if ~isfield(constrForceField{frame},'clusterAnalysis')
    display(['For frame: ',num2str(frame),' no cluster analysis has been performed!'])
    display('Nothing to do!')
    return;
end

global globYoung

globYoung=constrForceField{frame}.clusterAnalysis.par.globYoung;

u=constrForceField{frame}.clusterAnalysis.sol.u;
p=constrForceField{frame}.clusterAnalysis.mesh.p;
t=constrForceField{frame}.clusterAnalysis.mesh.t;

% output is interpolated to the node points p:
[u1,u2,exx,eyy,exy,sxx,syy,sxy,u1x,u2x,u1y,u2y]=postProcSol(u,p,t);

if nargin>4 && plotAll==1
    figure(3)
    pdesurf(p,t,u1)

    figure(4)
    pdesurf(p,t,u2)

    figure(5)
    pdesurf(p,t,sxx)

    figure(6)
    pdesurf(p,t,syy)

    figure(7)
    pdesurf(p,t,sxy)
end

currentImage = double(imread(imageFileList{frame}));
dPix_Plot=50;
min_x=min(constrForceField{frame}.roi.pos(:,1))-dPix_Plot;
max_x=max(constrForceField{frame}.roi.pos(:,1))+dPix_Plot;
min_y=min(constrForceField{frame}.roi.pos(:,2))-dPix_Plot;
max_y=max(constrForceField{frame}.roi.pos(:,2))+dPix_Plot;
[rows,cols]=size(currentImage);
j=1; % this arbitrary:
maxForcePlot=0.5*max(sqrt(constrForceField{frame}.clusterAnalysis.intf{j}.s_vec(:,1).^2+constrForceField{frame}.clusterAnalysis.intf{j}.s_vec(:,2).^2))/forceField(frame).par.gridSpacing;

%!!!% convert the force back to a stress to plot it with the
% traction stress. That is a bit dirty!
%factor_Pa2nN=constrForceField{frame}.par.gridSpacing^2*constrForceField{frame}.par.pixSize_mu^2*10^(-3); 
                
figure(8)
imagesc(constrForceField{frame}.clusterAnalysis.par.globYoung.val);
hold on
quiver(constrForceField{frame}.roi.pos(:,1),constrForceField{frame}.roi.pos(:,2),-constrForceField{frame}.roi.vec(:,1),-constrForceField{frame}.roi.vec(:,2),'-w');
plot(constrForceField{frame}.segmRes.curve(:,1),constrForceField{frame}.segmRes.curve(:,2),'k')
plot(constrForceField{frame}.segmRes.curveDilated(:,1),constrForceField{frame}.segmRes.curveDilated(:,2),'-b');
hold off
title(['Imposed Youngs modulus and inverted traction forces, frame: ',num2str(frame)])
axis equal
xlim([max([1 min_x]) min([cols,max_x])])
ylim([max([1 min_y]) min([rows,max_y])])
set(gca,'YDir','reverse')


figure(9)
imagesc(currentImage)
colormap('gray')
hold on
quiver(constrForceField{frame}.roi.pos(:,1),constrForceField{frame}.roi.pos(:,2),-constrForceField{frame}.roi.vec(:,1),-constrForceField{frame}.roi.vec(:,2),'-w');
plot(constrForceField{frame}.segmRes.curve(:,1),constrForceField{frame}.segmRes.curve(:,2),'k')
plot(constrForceField{frame}.segmRes.curveDilated(:,1),constrForceField{frame}.segmRes.curveDilated(:,2),'-b');
for j=1:length(constrForceField{frame}.network.edge)
    plot(constrForceField{frame}.network.edge{j}.intf(:,1),constrForceField{frame}.network.edge{j}.intf(:,2),'-k');
    plot(constrForceField{frame}.clusterAnalysis.intf{j}.cntrs(:,1),constrForceField{frame}.clusterAnalysis.intf{j}.cntrs(:,2),'-r');
    quiver(constrForceField{frame}.clusterAnalysis.intf{j}.cntrs(:,1),constrForceField{frame}.clusterAnalysis.intf{j}.cntrs(:,2),constrForceField{frame}.clusterAnalysis.intf{j}.s_vec(:,1)/maxForcePlot,constrForceField{frame}.clusterAnalysis.intf{j}.s_vec(:,2)/maxForcePlot,0)
    quiver(constrForceField{frame}.clusterAnalysis.bnd.pos(:,1),constrForceField{frame}.clusterAnalysis.bnd.pos(:,2),constrForceField{frame}.clusterAnalysis.bnd.s_vec(:,1)/maxForcePlot,constrForceField{frame}.clusterAnalysis.bnd.s_vec(:,2)/maxForcePlot,0)
end
hold off
title(['Stress profile along interface, frame: ',num2str(frame)])
axis equal
xlim([max([1 min_x]) min([cols,max_x])])
ylim([max([1 min_y]) min([rows,max_y])])
set(gca,'YDir','reverse')

figure(10)
imagesc(currentImage)
colormap('gray')
hold on
% quiver(constrForceField{frame}.roi.pos(:,1),constrForceField{frame}.roi.pos(:,2),-constrForceField{frame}.roi.vec(:,1),-constrForceField{frame}.roi.vec(:,2),'-w');
quiver(p(1,1:100:end),p(2,1:100:end),u1(1:100:end)',u2(1:100:end)','-m');
plot(constrForceField{frame}.segmRes.curve(:,1),constrForceField{frame}.segmRes.curve(:,2),'k')
plot(constrForceField{frame}.segmRes.curveDilated(:,1),constrForceField{frame}.segmRes.curveDilated(:,2),'-b');
for j=1:length(constrForceField{frame}.network.edge)
    plot(constrForceField{frame}.network.edge{j}.intf(:,1),constrForceField{frame}.network.edge{j}.intf(:,2),'-k');
    plot(constrForceField{frame}.clusterAnalysis.intf{j}.cntrs(:,1),constrForceField{frame}.clusterAnalysis.intf{j}.cntrs(:,2),'-r');
    quiver(constrForceField{frame}.clusterAnalysis.intf{j}.cntrs(:,1),constrForceField{frame}.clusterAnalysis.intf{j}.cntrs(:,2),constrForceField{frame}.clusterAnalysis.intf{j}.s_vec(:,1)/maxForcePlot,constrForceField{frame}.clusterAnalysis.intf{j}.s_vec(:,2)/maxForcePlot,0)
    quiver(constrForceField{frame}.clusterAnalysis.bnd.pos(:,1),constrForceField{frame}.clusterAnalysis.bnd.pos(:,2),constrForceField{frame}.clusterAnalysis.bnd.s_vec(:,1)/maxForcePlot,constrForceField{frame}.clusterAnalysis.bnd.s_vec(:,2)/maxForcePlot,0)
end
hold off
title(['Stress profile along interface, frame: ',num2str(frame)])
axis equal
xlim([max([1 min_x]) min([cols,max_x])])
ylim([max([1 min_y]) min([rows,max_y])])
set(gca,'YDir','reverse')
