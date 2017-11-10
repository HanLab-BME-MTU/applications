% To plot and save cell junction forces overlaid with cell or E-cad images
function plotClusterAnalysis_RNedit20121101(constrForceField,forceField,imageFileList,frameList,plotAll, doSave)

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

if nargin < 4 || isempty(frameList)
    frameList=input('Which frame(s) should be analyzed [default: all]: ');
    if isempty(frameList)
        for iframe=1:length(constrForceField);
            if isfield(constrForceField{iframe},'clusterAnalysis')
                frameList=[frameList,iframe];
            else
                display(['disregard frame: ',num2str(iframe),' because no cluster analysis was done'])
            end
        end
    end
end

if nargin<5 || isempty(plotAll)
    plotAll=0;
end

if nargin<6 || isempty(doSave)
    doSave = input('Want to save the cell-cell stress overlaid images? Yes=1; No =[0] ');
    if isempty(doSave)
        doSave = 0;
    else 
        doSave = 1;
    end
end

%% Determine the plot boundaries and pos of scale bars.
dPix=100;
dPixX=150;
dPixY=50;
textSpace=30;

% scaling factor of your choice to make the vectors in the quiver plots
% longer. If set to one, the longest vector just fits into the diagonal of
% the gridsize. Usually, a factor of 1. yields vectors that are too small.
% 1.5 seems to be a reasonable value. Of course, the length of the scale
% bars are scaled too.
scale_TF_vecs=1.5;
scale_intf_stress_vecs=1.5;

% set the values for the scale bars:
lengthScaleBar_mu=5;
fxScaleBar_Pa=20000;
fyScaleBar_Pa=0;
intfStressScaleBar_nN_per_um=10;

%% find the maximum stress values for setting the proper scales in the quiver plots: 
currentImage = double(imread(imageFileList{frameList(1)}));
[rows,cols]=size(currentImage);
xmin=Inf;
xmax=0;
ymin=Inf;
ymax=0;
forceScale=0;
maxIntfStressPlot=0;
for frame=frameList    
    currXmin=min(constrForceField{frame}.roi.pos(:,1))-dPix;
    currXmax=max(constrForceField{frame}.roi.pos(:,1))+dPix;
    currYmin=min(constrForceField{frame}.roi.pos(:,2))-dPix;
    currYmax=max(constrForceField{frame}.roi.pos(:,2))+dPix;
    if currXmin<xmin
        xmin=currXmin;
    end
    if currXmax>xmax
        xmax=currXmax;
    end
    if currYmin<ymin
        ymin=currYmin;
    end
    if currYmax>ymax
        ymax=currYmax;
    end
    
    % Do the scaling of the quiver plot of the stress by myself:
    maxForceComp=max(abs(constrForceField{frame}.roi.vec(:)));
    currForceScale=maxForceComp/constrForceField{frame}.par.gridSpacing;
    if currForceScale>forceScale
        forceScale=currForceScale;
    end
    
    % run through all omterfaces to find one maxIntfStressPlot for all
    % frames, so that the scale doesn't change from frame to frame! 
    for j=1:length(constrForceField{frame}.clusterAnalysis.intf)
        if ~isempty(constrForceField{frame}.clusterAnalysis.intf{j})
            maxIntfStressPlot=max(maxIntfStressPlot,max(sqrt(constrForceField{frame}.clusterAnalysis.intf{j}.s_vec(:,1).^2+constrForceField{frame}.clusterAnalysis.intf{j}.s_vec(:,2).^2))/forceField(frame).par.gridSpacing);
        end
    end
end
xLimVal=[max([1 xmin]) min([cols,xmax])];
yLimVal=[max([1 ymin]) min([rows,ymax])];

% apply the arbitrary scaling factors to make quiver vectors longer:
forceScale=forceScale/scale_TF_vecs;
maxIntfStressPlot=maxIntfStressPlot/scale_intf_stress_vecs;

%% make the plots:
for frame = frameList
    % first test if the cluster analysis has been performed on this frame. For
    % example, there is only a single cell:
    if ~isfield(constrForceField{frame},'clusterAnalysis')
        display(['For frame: ',num2str(frame),' no cluster analysis has been performed!'])
        display('Nothing to do!')
        return;
    end
    
    global globYoung    
    globYoung=constrForceField{frame}.clusterAnalysis.par.globYoung;
    
    % calculate the size of um in pixel for a scale bar:
    pixSize_mu=forceField(frame).par.pixSize_mu;
    lengthScaleBar_pix=lengthScaleBar_mu/pixSize_mu;
    
    u=constrForceField{frame}.clusterAnalysis.sol.u;
    p=constrForceField{frame}.clusterAnalysis.mesh.p;
    t=constrForceField{frame}.clusterAnalysis.mesh.t;
    
    % output is interpolated to the node points p:
    [u1,u2,exx,eyy,exy,sxx,syy,sxy,u1x,u2x,u1y,u2y]=postProcSol(u,p,t);
    
    if plotAll==1
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
        
    %!!!% convert the force back to a stress to plot it with the
    % traction stress. That is a bit dirty!
    %factor_Pa2nN=constrForceField{frame}.par.gridSpacing^2*constrForceField{frame}.par.pixSize_mu^2*10^(-3);
    
    % to convert the interfacial stresses from Pa*Pixel to nN/um:
    fac_PaPix_2_nN_per_um=constrForceField{frame}.par.pixSize_mu*10^(-3);
    % as the quiverplots are scaled by the maximum stress/force they are
    % dimensionless. Here, no conversion factor has to be applied. Only
    % the proper value for the scale bar needs to be calculated!
    
    figure(8)
    imagesc(constrForceField{frame}.clusterAnalysis.par.globYoung.val);
    hold on
    quiver(constrForceField{frame}.roi.pos(:,1),constrForceField{frame}.roi.pos(:,2),-constrForceField{frame}.roi.vec(:,1)/forceScale,-constrForceField{frame}.roi.vec(:,2)/forceScale,0,'-w');
    plot(constrForceField{frame}.segmRes.curve(:,1),constrForceField{frame}.segmRes.curve(:,2),'k')
    plot(constrForceField{frame}.segmRes.curveDilated(:,1),constrForceField{frame}.segmRes.curveDilated(:,2),'-b');
    % The scale bar um/pix:
    plot([xLimVal(2)-lengthScaleBar_pix-dPixX  xLimVal(2)-dPixX], [yLimVal(2)-dPixY yLimVal(2)-dPixY],'w','LineWidth',3)
    text( xLimVal(2)-lengthScaleBar_pix-dPixX, yLimVal(2)-dPixY-textSpace,[num2str(lengthScaleBar_mu),' \mum'],'HorizontalAlignment','left','color', 'w','FontSize',16)
    % The scale bar for the stresses:
    quiver(xLimVal(2)-lengthScaleBar_pix-dPixX, yLimVal(2)-2*dPixY,fxScaleBar_Pa/forceScale,fyScaleBar_Pa/forceScale,0,'w','LineWidth',2,'MaxHeadSize',5)
    text(  xLimVal(2)-lengthScaleBar_pix-dPixX, yLimVal(2)-2*dPixY-textSpace,[num2str(fxScaleBar_Pa/1000),' kPa'],'HorizontalAlignment','left','color', 'w','FontSize',16)
    axis equal
    xlim(xLimVal)
    ylim(yLimVal)
    hold off
    title(['Imposed Youngs modulus and inverted traction forces, frame: ',num2str(frame)])
    set(gca,'YDir','reverse')
    % if doSave == 1
    %   saveas(gcf,['ImposedYoungsModulusAndInvertedTractionForcesFrame',num2str(frame),'.tif'],'tif')
    %   saveas(gcf,['ImposedYoungsModulusAndInvertedTractionForcesFrame',num2str(frame),'.eps'],'psc2')
    % end
    
    figure(9)
    imagesc(currentImage)
    colormap('gray')
    hold on
    quiver(constrForceField{frame}.roi.pos(:,1),constrForceField{frame}.roi.pos(:,2),-constrForceField{frame}.roi.vec(:,1)/forceScale,-constrForceField{frame}.roi.vec(:,2)/forceScale,0,'-g','lineWidth',2);
    plot(constrForceField{frame}.segmRes.curve(:,1),constrForceField{frame}.segmRes.curve(:,2),'k')
    plot(constrForceField{frame}.segmRes.curveDilated(:,1),constrForceField{frame}.segmRes.curveDilated(:,2),'-b','lineWidth',3);
    for j=1:length(constrForceField{frame}.network.edge)
        plot(constrForceField{frame}.network.edge{j}.intf(:,1),constrForceField{frame}.network.edge{j}.intf(:,2),'-r','lineWidth',2);
        plot(constrForceField{frame}.clusterAnalysis.intf{j}.cntrs(:,1),constrForceField{frame}.clusterAnalysis.intf{j}.cntrs(:,2),'-r','lineWidth',2);
        quiver(constrForceField{frame}.clusterAnalysis.intf{j}.cntrs(:,1),constrForceField{frame}.clusterAnalysis.intf{j}.cntrs(:,2),constrForceField{frame}.clusterAnalysis.intf{j}.s_vec(:,1)/maxIntfStressPlot,constrForceField{frame}.clusterAnalysis.intf{j}.s_vec(:,2)/maxIntfStressPlot,0,'lineWidth',2)
        quiver(constrForceField{frame}.clusterAnalysis.bnd.pos(:,1),constrForceField{frame}.clusterAnalysis.bnd.pos(:,2),constrForceField{frame}.clusterAnalysis.bnd.s_vec(:,1)/maxIntfStressPlot,constrForceField{frame}.clusterAnalysis.bnd.s_vec(:,2)/maxIntfStressPlot,0,'b','lineWidth',2)
    end
    % The scale bar um/pix:
    plot([xLimVal(2)-lengthScaleBar_pix-dPixX  xLimVal(2)-dPixX], [yLimVal(2)-dPixY yLimVal(2)-dPixY],'w','LineWidth',3)
    text( xLimVal(2)-lengthScaleBar_pix-dPixX, yLimVal(2)-dPixY-textSpace,[num2str(lengthScaleBar_mu),' \mum'],'HorizontalAlignment','left','color', 'w','FontSize',16)
    % The scale bar for the stresses:
    quiver(xLimVal(2)-lengthScaleBar_pix-dPixX, yLimVal(2)-2*dPixY,fxScaleBar_Pa/forceScale,fyScaleBar_Pa/forceScale,0,'g','LineWidth',2,'MaxHeadSize',5)
    text(  xLimVal(2)-lengthScaleBar_pix-dPixX, yLimVal(2)-2*dPixY-textSpace,[num2str(fxScaleBar_Pa/1000),' kPa'],'HorizontalAlignment','left','color', 'g','FontSize',16)
    % The scale bar for the line stress/forces:
    quiver(xLimVal(2)-lengthScaleBar_pix-dPixX, yLimVal(2)-3*dPixY,intfStressScaleBar_nN_per_um/(maxIntfStressPlot*fac_PaPix_2_nN_per_um),0,0,'b','LineWidth',2,'MaxHeadSize',5)
    text(  xLimVal(2)-lengthScaleBar_pix-dPixX, yLimVal(2)-3*dPixY-textSpace,[num2str(intfStressScaleBar_nN_per_um),' nN/\mum'],'HorizontalAlignment','left','color', 'b','FontSize',16)
    hold off
    title(['Stress profile along interface, frame: ',num2str(frame)])
    axis equal
    xlim(xLimVal)
    ylim(yLimVal)
    set(gca,'YDir','reverse')
    if doSave == 1
        saveas(gcf,['StressProfileAlongInterfaceFrame',num2str(frame),'.tif'],'tif')
        saveas(gcf,['StressProfileAlongInterfaceFrame',num2str(frame),'.eps'],'psc2')
    end
    
    
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
        quiver(constrForceField{frame}.clusterAnalysis.intf{j}.cntrs(:,1),constrForceField{frame}.clusterAnalysis.intf{j}.cntrs(:,2),constrForceField{frame}.clusterAnalysis.intf{j}.s_vec(:,1)/maxIntfStressPlot,constrForceField{frame}.clusterAnalysis.intf{j}.s_vec(:,2)/maxIntfStressPlot,0)
        quiver(constrForceField{frame}.clusterAnalysis.bnd.pos(:,1),constrForceField{frame}.clusterAnalysis.bnd.pos(:,2),constrForceField{frame}.clusterAnalysis.bnd.s_vec(:,1)/maxIntfStressPlot,constrForceField{frame}.clusterAnalysis.bnd.s_vec(:,2)/maxIntfStressPlot,0)
    end
    % The scale bar um/pix:
    plot([xLimVal(2)-lengthScaleBar_pix-dPixX  xLimVal(2)-dPixX], [yLimVal(2)-dPixY yLimVal(2)-dPixY],'w','LineWidth',3)
    text( xLimVal(2)-lengthScaleBar_pix-dPixX, yLimVal(2)-dPixY-textSpace,[num2str(lengthScaleBar_mu),' \mum'],'HorizontalAlignment','left','color', 'w','FontSize',16)
    % The scale bar for the line stress/forces:
    quiver(xLimVal(2)-lengthScaleBar_pix-dPixX, yLimVal(2)-3*dPixY,intfStressScaleBar_nN_per_um/(maxIntfStressPlot*fac_PaPix_2_nN_per_um),0,0,'b','LineWidth',2,'MaxHeadSize',5)
    text(  xLimVal(2)-lengthScaleBar_pix-dPixX, yLimVal(2)-3*dPixY-textSpace,[num2str(intfStressScaleBar_nN_per_um),' nN/\mum'],'HorizontalAlignment','left','color', 'b','FontSize',16)
    hold off
    title(['Stress profile along interface, frame: ',num2str(frame)])
    axis equal
    xlim(xLimVal)
    ylim(yLimVal)
    set(gca,'YDir','reverse')
end