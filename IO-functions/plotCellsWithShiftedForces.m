function []=plotCellsWithShiftedForces(inputFileList,forceField,target_dir,frameList,displField)

% max of the color scale:
cMax=[];%3000;
cMaxZoom=cMax;
% If cMax is empty, the program will take the maximal stress value.

%read in Stack of bead images:
if nargin < 1 || isempty(inputFileList)
   [filename, pathname] = uigetfile({'*.tif';'*.jpg';'*.png';'*.*'}, ...
       'Select First Image of the Stack to be overlaied');
   
   if ~ischar(filename) || ~ischar(pathname)
       return;
   end
   
   inputFileList = getFileStackNames([pathname filesep filename]);
else
    isValid = 1;
    for i = 1:numel(inputFileList)
        isValid = isValid && exist(inputFileList{i}, 'file');
    end
    if ~isValid
        error('Invalid input files.');
    end
end

if nargin < 2 || isempty(forceField)
   [filename, pathname] = uigetfile({'*.mat';'*.*'}, ...
       'Select force field to be used as overlay');
       %the vector field:
       fileStruct=load([pathname filename]);
       forceField=fileStruct.forceField;
end

%get the target directory:
if nargin < 3 || isempty(target_dir)
    target_dir = uigetdir('','Please select target directory');
end


nImages = numel(inputFileList);
nVecFields = length(forceField);

if nargin < 4 || isempty(frameList)
    frameList=1:min(nImages,nVecFields);
end

needToCalcShift=false;
for i=frameList
    if ~isfield(forceField(i),'posShifted') || isempty(forceField(i).posShifted)
        needToCalcShift=true;
        break;
    end
end

if needToCalcShift && (nargin<5 || isempty(displField))
   [filename, pathname] = uigetfile({'*.mat';'*.*'}, ...
       'Select displField to be used for shifting the forceField');
       %the vector field:
       fileStruct=load([pathname filesep filename]);
       displField=fileStruct.displField;
end

if needToCalcShift
    [forceField]=shiftForceField(forceField,displField,frameList);
end

padZeros=floor(log10(length(nImages)))+1;

for i=1:min(nImages,nVecFields)
    imgInfo=imfinfo(inputFileList{1});
    maxX(i)=getfield(imgInfo,'Width');
    maxY(i)=getfield(imgInfo,'Height');
end
theXlim=max(maxX);
theYlim=max(maxY);

% Set the same force scale for all frames!
forceScale=0;
for i=frameList
    % Do the scaling of the quiver plot of the stress by myself:
    maxForceComp=max(abs(forceField(i).vec(:)));
    currForceScale=maxForceComp/forceField(i).par.gridSpacing;
    if currForceScale>forceScale
        forceScale=currForceScale;
    end
end
forceScale=1/3*forceScale;
displScale=2;

% calculate the length of the scale bars:
lengthScaleBar_mu=3;
uScaleBar_mu=lengthScaleBar_mu/displScale;
fxScaleBar_Pa=1000;
fyScaleBar_Pa=0;

display('Only the first frame is plotted!')
for i=frameList
    I = double(imread(inputFileList{i}));
    
    % extra spacing from the image edge in pixel:
    dPix=25;
    textSpace=12;
      
    % calculate the size of um in pixel for a scale bar:
    pixSize_mu=forceField(i).par.pixSize_mu;
    lengthScaleBar_pix=lengthScaleBar_mu/pixSize_mu;
    uScaleBar_pix=uScaleBar_mu/pixSize_mu;
    % This is the original force field:    
    figure(1)
    colormap('gray');
    imagesc(I)
    hold on
    if nargin>=5 && ~isempty(displField(i))
        quiver(displField(i).pos(:,1),displField(i).pos(:,2),displField(i).vec(:,1)*displScale,displField(i).vec(:,2)*displScale,0,'b');
    end    
    quiver(forceField(i).pos(:,1),forceField(i).pos(:,2),forceField(i).vec(:,1)/forceScale,forceField(i).vec(:,2)/forceScale,0,'r')
    % The scale bar um/pix:
    plot([theXlim-lengthScaleBar_pix-dPix theXlim-dPix], [theYlim-dPix theYlim-dPix],'w','LineWidth',3)
    text(theXlim-lengthScaleBar_pix-dPix, theYlim-dPix-textSpace,[num2str(lengthScaleBar_mu),' \mum'],'HorizontalAlignment','left','color', 'w','FontSize',16)
    % The scale bar for the stresses:
    quiver(theXlim-lengthScaleBar_pix-dPix,theYlim-2*dPix,fxScaleBar_Pa/forceScale,fyScaleBar_Pa/forceScale,0,'w','LineWidth',2,'MaxHeadSize',5)
    text(theXlim-lengthScaleBar_pix-dPix, theYlim-2*dPix-textSpace,[num2str(fxScaleBar_Pa/1000),' kPa'],'HorizontalAlignment','left','color', 'w','FontSize',16)
    if nargin>=5 && ~isempty(displField(i))
        % The scale bar for the displacement:
        quiver(theXlim-lengthScaleBar_pix-dPix,theYlim-3*dPix,uScaleBar_pix*displScale,0,0,'w','LineWidth',2,'MaxHeadSize',5)
        text(theXlim-lengthScaleBar_pix-dPix, theYlim-3*dPix-textSpace,[num2str(uScaleBar_mu),' \mum'],'HorizontalAlignment','left','color', 'w','FontSize',16)
    end
    %!!! Equal axis is important for dimension/y-forces scale to be accuarte!!!
    axis equal
    xlim([1 theXlim])
    ylim([1 theYlim])
    title('Cells with original forces')
    set(gca,'YDir','reverse')%,'XTick',[],'YTick',[])
    %saveas(gcf,[target_dir,filesep,'Cells_with_',fieldName,num2str(i,['%0.',int2str(padZeros),'d']),'.tiff'],'tiffn');
    %saveas(gcf,[target_dir,filesep,'Cells_with_',fieldName,num2str(i,['%0.',int2str(padZeros),'d']),'.eps'], 'psc2');
    hold off
    
%     theXlim=457;
%     theYlim=457;
%     xrange = [45 theXlim];
%     yrange = [45 theYlim];
    
    xrange = [0 theXlim];
    yrange = [0 theYlim];
    
    % This is the shifted force field:
    figure(2)
    colormap('gray');
    imagesc(I)
    hold on
    if nargin>=5 && ~isempty(displField(i))
        quiver(displField(i).pos(:,1),displField(i).pos(:,2),displField(i).vec(:,1)*displScale,displField(i).vec(:,2)*displScale,0,'b');
    end
    quiver(forceField(i).posShifted(:,1),forceField(i).posShifted(:,2),forceField(i).vec(:,1)/forceScale,forceField(i).vec(:,2)/forceScale,0,'r')
    % The scale bar um/pix:
    plot([theXlim-lengthScaleBar_pix-dPix theXlim-dPix], [theYlim-dPix theYlim-dPix],'w','LineWidth',3)
    text(theXlim-lengthScaleBar_pix-dPix, theYlim-dPix-textSpace,[num2str(lengthScaleBar_mu),' \mum'],'HorizontalAlignment','left','color', 'w','FontSize',16)
    % The scale bar for the stresses:
    quiver(theXlim-lengthScaleBar_pix-dPix,theYlim-2*dPix,fxScaleBar_Pa/forceScale,fyScaleBar_Pa/forceScale,0,'w','LineWidth',2,'MaxHeadSize',5)
    text(theXlim-lengthScaleBar_pix-dPix, theYlim-2*dPix-textSpace,[num2str(fxScaleBar_Pa/1000),' kPa'],'HorizontalAlignment','left','color', 'w','FontSize',16)
    if nargin>=5 && ~isempty(displField(i))
        % The scale bar for the displacement:
        quiver(theXlim-lengthScaleBar_pix-dPix,theYlim-3*dPix,uScaleBar_pix*displScale,0,0,'w','LineWidth',2,'MaxHeadSize',5)
        text(theXlim-lengthScaleBar_pix-dPix, theYlim-3*dPix-textSpace,[num2str(uScaleBar_mu),' \mum'],'HorizontalAlignment','left','color', 'w','FontSize',16)
    end
    %!!! Equal axis is important for dimension/y-forces scale to be accuarte!!!
    axis equal
    title(['Cells with shifted forces (red) and displ. (blue), reg. param=',num2str(forceField(i).par.regParam)])
    set(gca,'YDir','reverse') %,'XTick',[],'YTick',[])
    %xlim([1 theXlim])
    %ylim([1 theYlim])
    xlim(xrange)
    ylim(yrange)
    saveas(gcf,[target_dir,filesep,'Cells_with_shifted_forces',num2str(i,['%0.',int2str(padZeros),'d']),'.tiff'],'tiffn');
    saveas(gcf,[target_dir,filesep,'Cells_with_shifted_forces',num2str(i,['%0.',int2str(padZeros),'d']),'.eps'], 'psc2');
    % print('-depsc2','-loose', [target_dir,filesep,'Cells_with_shifted_forces_HQ',num2str(i,['%0.',int2str(padZeros),'d']),'.eps']);
    hold off
    
    
    % This is the direct comparison of the two force fields:
    figure(3)
    colormap('gray');
    imagesc(I)
    hold on
    quiver(forceField(i).pos(:,1),forceField(i).pos(:,2),forceField(i).vec(:,1)/forceScale,forceField(i).vec(:,2)/forceScale,0,'b')
    quiver(forceField(i).posShifted(:,1),forceField(i).posShifted(:,2),forceField(i).vec(:,1)/forceScale,forceField(i).vec(:,2)/forceScale,0,'r')
    % The scale bar um/pix:
    plot([theXlim-lengthScaleBar_pix-dPix theXlim-dPix], [theYlim-dPix theYlim-dPix],'w','LineWidth',3)
    text(theXlim-lengthScaleBar_pix-dPix, theYlim-dPix-textSpace,[num2str(lengthScaleBar_mu),' \mum'],'HorizontalAlignment','left','color', 'w','FontSize',16)
    % The scale bar for the stresses:
    quiver(theXlim-lengthScaleBar_pix-dPix,theYlim-2*dPix,fxScaleBar_Pa/forceScale,fyScaleBar_Pa/forceScale,0,'w','LineWidth',2,'MaxHeadSize',5)
    text(theXlim-lengthScaleBar_pix-dPix, theYlim-2*dPix-textSpace,[num2str(fxScaleBar_Pa/1000),' kPa'],'HorizontalAlignment','left','color', 'w','FontSize',16)
%!!! Equal axis is important for dimension/y-forces scale to be accuarte!!!
    axis equal 
    xlim(xrange)
    ylim(yrange)
    title('Comparison of the two force fields')
    set(gca,'YDir','reverse')%,'XTick',[],'YTick',[])
    %saveas(gcf,[target_dir,filesep,'Cells_with_',fieldName,num2str(i,['%0.',int2str(padZeros),'d']),'.tiff'],'tiffn');
    %saveas(gcf,[target_dir,filesep,'Cells_with_',fieldName,num2str(i,['%0.',int2str(padZeros),'d']),'.eps'], 'psc2');
    hold off
    
    % cut off for evaluating background:
    xLimZoom=input('Input xLim for zoom in [300 430]:');
    if isempty(xLimZoom)
        xLimZoom=[300 430];
    end
    
    yLimZoom=input('Input yLim for zoom in [100 230]:');
    if isempty(yLimZoom)
        yLimZoom=[100 230];
    end
    
    % cut off for evaluating background:
    xCutOff=input('Input cutoff for background evaluation [450]:');
    if isempty(xCutOff)
        xCutOff=450;
    end
    
    % Create a Mask from the Zoom in region:
    bwMask(yLimZoom(1):yLimZoom(2),xLimZoom(1):xLimZoom(2))=1;
    [inpos,invec]=findVectorFieldInMask(forceField(i).pos,forceField(i).vec,bwMask);
    forceFieldZoom(i).pos=inpos;
    forceFieldZoom(i).vec=invec;   
    
    [forceFieldZoom]=shiftForceField(forceFieldZoom,displField,frameList); 
    
    
    % Do the scaling of the quiver plot of the stress by myself:
    maxForceCompZoom=max(abs(forceFieldZoom(i).vec(:)));
    forceScaleZoom=maxForceCompZoom/forceField(i).par.gridSpacing;
    
    % Spacings for the scale bars:
    dPixZoom=20;
    textSpaceZoom=5;
    
    % calculate the size of um in pixel for a scale bar:
    lengthScaleBarZoom_mu=2;
    lengthScaleBarZoom_pix=lengthScaleBarZoom_mu/pixSize_mu;
   
    % calculate the scale bar for the stress vectors:
    fxScaleBarZoom_Pa=2000;
    fyScaleBarZoom_Pa=0;
    
    figure(4)
    colormap('gray');
    imagesc(I)
    hold on
    quiver(forceFieldZoom(i).posShifted(:,1),forceFieldZoom(i).posShifted(:,2),forceFieldZoom(i).vec(:,1)/forceScaleZoom,forceFieldZoom(i).vec(:,2)/forceScaleZoom,0,'y')
    % The scale bar um/pix:
    plot([xLimZoom(2)-lengthScaleBarZoom_pix-dPixZoom xLimZoom(2)-dPixZoom], [yLimZoom(2)-dPixZoom yLimZoom(2)-dPixZoom],'w','LineWidth',3)
    text(xLimZoom(2)-lengthScaleBarZoom_pix-dPixZoom, yLimZoom(2)-dPixZoom-textSpaceZoom,[num2str(lengthScaleBarZoom_mu),' \mum'],'HorizontalAlignment','left','color', 'w','FontSize',16)
    % The scale bar for the stresses:
    quiver(xLimZoom(2)-lengthScaleBarZoom_pix-dPixZoom,yLimZoom(2)-2*dPixZoom,fxScaleBarZoom_Pa/forceScaleZoom,fyScaleBarZoom_Pa/forceScaleZoom,0,'w','LineWidth',2,'MaxHeadSize',5)
    text(xLimZoom(2)-lengthScaleBarZoom_pix-dPixZoom, yLimZoom(2)-2*dPixZoom-textSpaceZoom,[num2str(fxScaleBarZoom_Pa/1000),' kPa'],'HorizontalAlignment','left','color', 'w','FontSize',16)
%!!! Equal axis is important for dimension/y-forces scale to be accuarte!!!
    axis equal 
    xlim(xLimZoom)
    ylim(yLimZoom)
    title('Force field shifted and zoomed')
    set(gca,'YDir','reverse')%,'XTick',[],'YTick',[])
    %saveas(gcf,[target_dir,filesep,'Cells_with_',fieldName,num2str(i,['%0.',int2str(padZeros),'d']),'.tiff'],'tiffn');
    %saveas(gcf,[target_dir,filesep,'Cells_with_',fieldName,num2str(i,['%0.',int2str(padZeros),'d']),'.eps'], 'psc2');
    hold off    
    
    
    maxStress=max(sqrt(sum(forceField(i).vec.^2,2)));
    if isempty(cMax);
        cMax=maxStress;
    end
  
    [rows,cols]=size(I);
    [XI,YI]=meshgrid(1:cols,1:rows);
    % for the magnitude:
    Mblue = griddata(forceField(i).posShifted(:,1),forceField(i).posShifted(:,2),sqrt(sum(forceField(i).vec.^2,2)),XI,YI,'cubic');
    % for x-component:
    % Mblue = griddata(forceField(i).posShifted(:,1),forceField(i).posShifted(:,2),forceField(i).vec(:,1),XI,YI,'cubic');
    % for y-component:
    % Mblue = griddata(forceField(i).posShifted(:,1),forceField(i).posShifted(:,2),forceField(i).vec(:,2),XI,YI,'cubic');
    % remove NaNs:
    Mblue(isnan(Mblue))=0;
    % plot colorful image:
    figure(5)
    colormap('jet')
    imagesc(Mblue)
    hold on
    % The scale bar um/pix:
    plot([theXlim-lengthScaleBar_pix-dPix theXlim-dPix], [theYlim-dPix theYlim-dPix],'w','LineWidth',3)
    text(theXlim-lengthScaleBar_pix-dPix, theYlim-dPix-textSpace,[num2str(lengthScaleBar_mu),' \mum'],'HorizontalAlignment','left','color', 'w','FontSize',16)
    title('Force magnitude raw')
    caxis([0 cMax])
%!!! Equal axis is important for dimension/y-forces scale to be accuarte!!!
    axis equal    
    xlim(xrange)
    ylim(yrange)
    hold off
    
    
    maxStressZoom=max(sqrt(sum(forceFieldZoom(i).vec.^2,2)));
    if isempty(cMaxZoom);
        cMaxZoom=maxStressZoom;
    end    
    
    figure(6)
    colormap('jet')
    imagesc(Mblue)
    hold on
    % The scale bar um/pix:
    plot([xLimZoom(2)-lengthScaleBarZoom_pix-dPixZoom xLimZoom(2)-dPixZoom], [yLimZoom(2)-dPixZoom yLimZoom(2)-dPixZoom],'w','LineWidth',3)
    text(xLimZoom(2)-lengthScaleBarZoom_pix-dPixZoom, yLimZoom(2)-dPixZoom-textSpaceZoom,[num2str(lengthScaleBarZoom_mu),' \mum'],'HorizontalAlignment','left','color', 'w','FontSize',16)
%!!! Equal axis is important for dimension/y-forces scale to be accuarte!!!
    axis equal 
    xlim(xLimZoom)
    ylim(yLimZoom)
    title('Force magnitude raw zoomed')
    caxis([0 cMaxZoom])
    %saveas(gcf,[target_dir,filesep,'Cells_with_',fieldName,num2str(i,['%0.',int2str(padZeros),'d']),'.tiff'],'tiffn');
    %saveas(gcf,[target_dir,filesep,'Cells_with_',fieldName,num2str(i,['%0.',int2str(padZeros),'d']),'.eps'], 'psc2');
    hold off
    
    
    % magnitude and vectors:
    cMaxDispl=max(sqrt(sum(displField(i).vec.^2,2)));
    MblueDispl = griddata(displField(i).pos(:,1),displField(i).pos(:,2),sqrt(sum(displField(i).vec.^2,2)),XI,YI,'cubic');
    % remove NaNs:
    MblueDispl(isnan(MblueDispl))=0;
    
    figure(7)
    colormap('jet')
    imagesc(MblueDispl*pixSize_mu)
    hold on
    if nargin>=5 && ~isempty(displField(i))
        quiver(displField(i).pos(:,1),displField(i).pos(:,2),displField(i).vec(:,1)*displScale,displField(i).vec(:,2)*displScale,0,'w');
    end
    % The scale bar um/pix:
    plot([theXlim-lengthScaleBar_pix-dPix theXlim-dPix], [theYlim-dPix theYlim-dPix],'w','LineWidth',3)
    text(theXlim-lengthScaleBar_pix-dPix, theYlim-dPix-textSpace,[num2str(lengthScaleBar_mu),' \mum'],'HorizontalAlignment','left','color', 'w','FontSize',16)
    if nargin>=5 && ~isempty(displField(i))
        % The scale bar for the displacement:
        quiver(theXlim-lengthScaleBar_pix-dPix,theYlim-3*dPix,uScaleBar_pix*displScale,0,0,'w','LineWidth',2,'MaxHeadSize',5)
        text(theXlim-lengthScaleBar_pix-dPix, theYlim-3*dPix-textSpace,[num2str(uScaleBar_mu),' \mum'],'HorizontalAlignment','left','color', 'w','FontSize',16)
    end
    %!!! Equal axis is important for dimension/y-forces scale to be accuarte!!!
    axis equal
    set(gca,'YDir','reverse') %,'XTick',[],'YTick',[])
    %xlim([1 theXlim])
    %ylim([1 theYlim])
    title('Displ magnitude and vec')
    xlim(xrange)
    ylim(yrange)
    hold off
    
    figure(8)
    colormap('jet')
    imagesc(Mblue)
    hold on
    % The scale bar um/pix:
    plot([theXlim-lengthScaleBar_pix-dPix theXlim-dPix], [theYlim-dPix theYlim-dPix],'w','LineWidth',3)
    text(theXlim-lengthScaleBar_pix-dPix, theYlim-dPix-textSpace,[num2str(lengthScaleBar_mu),' \mum'],'HorizontalAlignment','left','color', 'w','FontSize',16)
    caxis([0 cMax])
    quiver(forceField(i).posShifted(:,1),forceField(i).posShifted(:,2),forceField(i).vec(:,1)/forceScale,forceField(i).vec(:,2)/forceScale,0,'w')
    % The scale bar um/pix:
    plot([theXlim-lengthScaleBar_pix-dPix theXlim-dPix], [theYlim-dPix theYlim-dPix],'w','LineWidth',3)
    text(theXlim-lengthScaleBar_pix-dPix, theYlim-dPix-textSpace,[num2str(lengthScaleBar_mu),' \mum'],'HorizontalAlignment','left','color', 'w','FontSize',16)
    % The scale bar for the stresses:
    quiver(theXlim-lengthScaleBar_pix-dPix,theYlim-2*dPix,fxScaleBar_Pa/forceScale,fyScaleBar_Pa/forceScale,0,'w','LineWidth',2,'MaxHeadSize',5)
    text(theXlim-lengthScaleBar_pix-dPix, theYlim-2*dPix-textSpace,[num2str(fxScaleBar_Pa/1000),' kPa'],'HorizontalAlignment','left','color', 'w','FontSize',16)
    %!!! Equal axis is important for dimension/y-forces scale to be accuarte!!!
    axis equal
    set(gca,'YDir','reverse') %,'XTick',[],'YTick',[])
    %xlim([1 theXlim])
    %ylim([1 theYlim])
    title(['Force magnitude and vec (raw) reg. param=',num2str(forceField(i).par.regParam)])
    xlim(xrange)
    ylim(yrange)
    hold off
    
    figure(9)
    colormap('gray');
    imagesc(max(Mblue(:))*I/max(I(:)))
    hold on
    %contour(max(I(:))*Mblue/max(Mblue(:)),10);
    contour(Mblue,10);
    % The scale bar um/pix:
    plot([theXlim-lengthScaleBar_pix-dPix theXlim-dPix], [theYlim-dPix theYlim-dPix],'w','LineWidth',3)
    text(theXlim-lengthScaleBar_pix-dPix, theYlim-dPix-textSpace,[num2str(lengthScaleBar_mu),' \mum'],'HorizontalAlignment','left','color', 'w','FontSize',16)
    xlim([1 theXlim])
    ylim([1 theYlim])
    title('Force magnitude raw')
    set(gca,'YDir','reverse')%,'XTick',[],'YTick',[])
%!!! Equal axis is important for dimension/y-forces scale to be accuarte!!!
    axis equal    
    xlim(xrange)
    ylim(yrange)
    hold off
    
    % substract background:
    checkVec=forceField(i).posShifted(:,1)<xCutOff;
    backGround(i)=forceField(i);
    backGround(i).pos(checkVec,:)=[];
    backGround(i).vec(checkVec,:)=[];
    bgLevel=mean(sqrt(sum(backGround(i).vec.^2,2)));    
    Mblue=Mblue-bgLevel;
    Mblue(Mblue<0)=0;
    Mblue  =Mblue/max(Mblue(:));
    
    figure(10)
    colormap('jet')
    imagesc(Mblue)    
    hold on;
    caxis([0 1])
    title('Force magnitude background substracted')
    set(gca,'YDir','reverse')%,'XTick',[],'YTick',[])    
    % The scale bar um/pix:
    plot([theXlim-lengthScaleBar_pix-dPix theXlim-dPix], [theYlim-dPix theYlim-dPix],'w','LineWidth',3)
    text(theXlim-lengthScaleBar_pix-dPix, theYlim-dPix-textSpace,[num2str(lengthScaleBar_mu),' \mum'],'HorizontalAlignment','left','color', 'w','FontSize',16)
%!!! Equal axis is important for dimension/y-forces scale to be accuarte!!!
    axis equal    
    xlim(xrange)
    ylim(yrange)
    
    
    max(Mblue(:))
    
    
    figure(11)
    colormap('gray');
    imagesc(I)
    hold on
    maxI=max(I(:));
    contour(maxI*Mblue/max(Mblue(:)),13);
    % The scale bar um/pix:
    plot([theXlim-lengthScaleBar_pix-dPix theXlim-dPix], [theYlim-dPix theYlim-dPix],'w','LineWidth',3)
    text(theXlim-lengthScaleBar_pix-dPix, theYlim-dPix-textSpace,[num2str(lengthScaleBar_mu),' \mum'],'HorizontalAlignment','left','color', 'w','FontSize',16)
    hold off
    xlim([1 theXlim])
    ylim([1 theYlim])
    title('Force magnitude background substracted')
    set(gca,'YDir','reverse')%,'XTick',[],'YTick',[])
%!!! Equal axis is important for dimension/y-forces scale to be accuarte!!!
    axis equal    
    xlim(xrange)
    ylim(yrange)
    
    
    Mred  =I/max(I(:));
    Mgreen=zeros(size(I));
    %Mblue =zeros(size(I));
    
    Mblue  =Mblue/max(Mblue(:));    
    Mblue(Mblue(:)<0)=0;
    
    RGBmat(1:rows,1:cols,1)=abs(Mred-1);
    RGBmat(1:rows,1:cols,2)=Mgreen;
    RGBmat(1:rows,1:cols,3)=Mblue;
    
    figure(12)
    imagesc(RGBmat)
    hold on
    % The scale bar um/pix:
    plot([theXlim-lengthScaleBar_pix-dPix theXlim-dPix], [theYlim-dPix theYlim-dPix],'w','LineWidth',3)
    text(theXlim-lengthScaleBar_pix-dPix, theYlim-dPix-textSpace,[num2str(lengthScaleBar_mu),' \mum'],'HorizontalAlignment','left','color', 'w','FontSize',16)
    title('Cells with force Magnitude')
%!!! Equal axis is important for dimension/y-forces scale to be accuarte!!!
    axis equal    
    xlim(xrange)
    ylim(yrange)
    %imwrite(RGBmat,[target_dir,filesep,'Cells_with_',fieldName,'Magnitude',num2str(i,['%0.',int2str(padZeros),'d']),'.tiff']) 
    hold off
end