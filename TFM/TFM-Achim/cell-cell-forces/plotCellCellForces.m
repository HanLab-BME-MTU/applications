function plotCellCellForces(constrForceField,forceField,trackedNet,toDoList,imageFileList,target_dir,showTFHeatMap,showIntfStress,doSave,movieFormat,fps)
% This function plots the residual forces in the cell cluster!
% Cells are tracked (maintain the same color code) only if trackedNet is 
% provided. In case trackedNet is not provided the color code will 
% most like not be consistent throughout the movie.
%
% last change: 03/24/2013, implementing constant color code for cells.
% the same could be done for the interfaces to (with the same trick) but
% not yet needed.

if nargin < 1 || isempty(constrForceField)
    [filename_constrForceField, pathname_constrForceField] = uigetfile({'*.mat';'*.*'}, ...
        'Select cellCellForces.mat to be used');
    fileStruct=load([pathname_constrForceField filesep filename_constrForceField]);
    constrForceField=fileStruct.constrForceField;
    
    target_dir=pathname_constrForceField;
end

if nargin < 2 || isempty(forceField)
    [filename_forceField, pathname_forceField] = uigetfile({'*.mat';'*.*'}, ...
        'Select forceField.mat to be used');
    fileStruct=load([pathname_forceField filesep filename_forceField]);
    forceField=fileStruct.forceField;
end


if nargin < 3 || isempty(trackedNet)
    try
        [filename_trackedNet, pathname_trackedNet] = uigetfile({'*.mat';'*.*'}, ...
            'Select trackedNet.mat to be used');
        fileStruct=load([pathname_trackedNet filesep filename_trackedNet]);
        trackedNet=fileStruct.trackedNet;
    catch
        trackedNet=[];
    end
end

%  ------
if nargin<4 || isempty(toDoList)
    toDoList=[];
    for frame=1:length(constrForceField)
        if ~isempty(constrForceField{frame}) && isfield(constrForceField{frame},'cell')
            toDoList=horzcat(toDoList,frame);
        end
    end
end

%read in Stack of images:
if nargin < 5 || isempty(imageFileList)
    [filename, pathname] = uigetfile({'*.tif';'*.jpg';'*.png';'*.*'}, ...
        'Select First Ecad-Image');
    
    if ~ischar(filename) || ~ischar(pathname)
        return;
    end
    
    imageFileList = getFileStackNames([pathname filesep filename]);
else
    isValid = 1;
    for frame = 1:numel(imageFileList)
        isValid = isValid && exist(imageFileList{frame}, 'file');
    end
    if ~isValid
        error('Invalid input files.');
    end
end

%get the target directory:
if nargin < 6 || isempty(target_dir)
    target_dir = uigetdir('','Please select target directory');
end

% the default is showing the Ecad-Images:
if nargin < 7 || isempty(showTFHeatMap)
    showTFHeatMap=0;
end

% the default is that the interface-stresses are shown:
if nargin < 8 || isempty(showIntfStress)
    showIntfStress=1;
end

if nargin < 9 || isempty(doSave)
    doSave = 1;
end

if nargin < 10 || isempty(fps)
    fps = 8;
end

movieFileName='mov_CellCellForces';
if nargin == 10 && (isempty(movieFormat) || strcmp(movieFormat,'mov') == 1 || strcmp(movieFormat,'MOV'))
    movieFileName = [movieFileName,'.mov'];
    MakeQTMovie('start',movieFileName);
    MakeQTMovie('framerate',fps);
    MakeQTMovie('quality',1);
    movieFormat ='mov';
    doMovie=1;
elseif nargin == 10 && (strcmp(movieFormat,'avi') == 1 || strcmp(movieFormat,'AVI'))
    movieFileName = [movieFileName,'.avi'];
    movieFormat ='avi';
    doMovie=1;
else
    doMovie=0;
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
currentImage = double(imread(imageFileList{toDoList(1)}));
[rows,cols]=size(currentImage);
xmin=Inf;
xmax=0;
ymin=Inf;
ymax=0;
forceScale=0;
maxIntfStressPlot=0;
cMax=-Inf;
for frame=toDoList
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
    
    % only needed for the TFM-heat-map:
    if showTFHeatMap
        currMaxForce=max(sqrt(sum(constrForceField{frame}.roi.vec.^2,2)));
        cMax = max(cMax,currMaxForce);
    end
    
    % only needed for showing the stresses along the interfaces:
    if showIntfStress
        % run through all interfaces to find one maxIntfStressPlot for all
        % frames, so that the scale doesn't change from frame to frame! 
        for j=1:length(constrForceField{frame}.clusterAnalysis.intf)
            if ~isempty(constrForceField{frame}.clusterAnalysis.intf{j})
                maxIntfStressPlot=max(maxIntfStressPlot,max(sqrt(constrForceField{frame}.clusterAnalysis.intf{j}.s_vec(:,1).^2+constrForceField{frame}.clusterAnalysis.intf{j}.s_vec(:,2).^2))/forceField(frame).par.gridSpacing);
            end
        end
    end
end
%% apply the arbitrary scaling factors to make quiver vectors longer:
forceScale=forceScale/scale_TF_vecs;
maxIntfStressPlot=maxIntfStressPlot/scale_intf_stress_vecs;

xLimVal=[max([1 xmin]) min([cols,xmax])];
yLimVal=[max([1 ymin]) min([rows,ymax])];


%% plot:
for frame=toDoList
    display(['work on frame: ',num2str(frame)]);
    currentImage = double(imread(imageFileList{frame}));
    
    % calculate the size of um in pixel for a scale bar:
    pixSize_mu=forceField(frame).par.pixSize_mu;
    lengthScaleBar_pix=lengthScaleBar_mu/pixSize_mu;
    
    
    % Number of cells and interfaces to plot:
    numCells=length(constrForceField{frame}.cell);
    numIntf =length(constrForceField{frame}.interface);
    
    % scaling factor to convert residual forces [nN] back to stresses:
    factor_Pa2nN=constrForceField{frame}.par.gridSpacing^2*constrForceField{frame}.par.pixSize_mu^2*10^(-3);
    
    % to convert the interfacial stresses from Pa*Pixel to nN/um:
    fac_PaPix_2_nN_per_um=constrForceField{frame}.par.pixSize_mu*10^(-3);
    
    % Here starts the image:
    scrsz = get(0,'ScreenSize');
    h=figure('Position',[1 scrsz(4) scrsz(3) scrsz(4)]);
    if showTFHeatMap
        [XI,YI]=meshgrid(1:cols,1:rows);
        % for plotting the TFM magnitude at the unstrained position:
        Mblue  = griddata(forceField(frame).pos(:,1),forceField(frame).pos(:,2),sqrt(sum(forceField(frame).vec.^2,2)),XI,YI,'cubic');
        % for plotting the TFM magnitude at the deformed position:
        % Mblue  = griddata(forceField(frame).posShifted(:,1),forceField(frame).posShifted(:,2),sqrt(sum(forceField(frame).vec.^2,2)),XI,YI,'cubic');
        % remove NaNs:
        Mblue(isnan(Mblue))=0;
        colormap('jet')
        imagesc(Mblue); 
        caxis([0 cMax]); cb=colorbar; ylab=get(cb,'ylabel');
        set(ylab,'String','traction stress magnitude [Pa]'); 
    else
        imagesc(currentImage)
        colormap('gray')
    end    
    hold on
    text((xLimVal(2)-xLimVal(1))*0.05,(yLimVal(2)-yLimVal(1))*0.05,num2str(frame),'Color','white','FontSize',18);
    %plot inner boundary in black:
    plot(constrForceField{frame}.segmRes.curve(:,1),constrForceField{frame}.segmRes.curve(:,2),'k')
    %plot complete forceField:
    %quiver(forceField(frame).pos(:,1),forceField(frame).pos(:,2),forceField(frame).vec(:,1)/forceScale,forceField(frame).vec(:,2)/forceScale,0,'g');
    marker=['r','b','m','c','g','y','w'];
    % plot the traction force field:
    for k=1:numCells
        % find tracked cellID:
        if length(trackedNet)<frame || isempty(trackedNet{frame})
            display('Tracking not possible: trackedNet is empty')
            cellIdTracked=k;
        else
            cellIdTracked=findCellIdTrackedNet(trackedNet{frame},constrForceField{frame}.cell{k}.center);
        end
        %plot cell{j}:
        quiver(constrForceField{frame}.cell{k}.pos(:,1),constrForceField{frame}.cell{k}.pos(:,2),constrForceField{frame}.cell{k}.vec(:,1)/forceScale,constrForceField{frame}.cell{k}.vec(:,2)/forceScale,0,marker(mod(cellIdTracked,7)+1));
        %!!!    % convert the force back to a stress to plot it with the
        % traction stress. That is a bit dirty!
        resForcePos=constrForceField{frame}.cell{k}.stats.resForce.pos;
        quiver(resForcePos(1),resForcePos(2),constrForceField{frame}.cell{k}.stats.resForce.vec(:,1)/(forceScale*factor_Pa2nN),constrForceField{frame}.cell{k}.stats.resForce.vec(:,2)/(forceScale*factor_Pa2nN),0,marker(mod(cellIdTracked,7)+1));
        plot(resForcePos(1),resForcePos(2),['o',marker(mod(cellIdTracked,7)+1)])
        %plot boundaries:
        plot(constrForceField{frame}.cell{k}.boundary(:,1),constrForceField{frame}.cell{k}.boundary(:,2),marker(mod(cellIdTracked,7)+1))
    end
    % plot the pixel-precise interface:
    for k=1:numIntf
        plot(constrForceField{frame}.interface{k}.pos(:,1),constrForceField{frame}.interface{k}.pos(:,2),'--w')
    end
    % plot the holes:
    if isfield(constrForceField{frame}.segmRes,'hole')
        for holeId=1:length(constrForceField{frame}.segmRes.hole)
            plot(constrForceField{frame}.segmRes.hole{holeId}.curve(:,1),constrForceField{frame}.segmRes.hole{holeId}.curve(:,2),'k')
        end
    end
    % plot the stress along the interfaces if wanted:
    if showIntfStress
        for j=1:length(constrForceField{frame}.network.edge)
            % for highlighting the sparse interface (centers are only every tenth pixels) with a white line:
            plot(constrForceField{frame}.clusterAnalysis.intf{j}.cntrs(:,1),constrForceField{frame}.clusterAnalysis.intf{j}.cntrs(:,2),'-w','lineWidth',2);
            quiver(constrForceField{frame}.clusterAnalysis.intf{j}.cntrs(:,1),constrForceField{frame}.clusterAnalysis.intf{j}.cntrs(:,2),constrForceField{frame}.clusterAnalysis.intf{j}.s_vec(:,1)/maxIntfStressPlot,constrForceField{frame}.clusterAnalysis.intf{j}.s_vec(:,2)/maxIntfStressPlot,0,'w','lineWidth',2)
        end
        % The scale bar for the line stress/forces:
        quiver(xLimVal(2)-lengthScaleBar_pix-dPixX, yLimVal(2)-3*dPixY,intfStressScaleBar_nN_per_um/(maxIntfStressPlot*fac_PaPix_2_nN_per_um),0,0,'w','LineWidth',2,'MaxHeadSize',5)
        text(  xLimVal(2)-lengthScaleBar_pix-dPixX, yLimVal(2)-3*dPixY-textSpace,[num2str(intfStressScaleBar_nN_per_um),' nN/\mum'],'HorizontalAlignment','left','color', 'w','FontSize',16)
    end
    % The scale bar um/pix:
    plot([xLimVal(2)-lengthScaleBar_pix-dPixX  xLimVal(2)-dPixX], [yLimVal(2)-dPixY yLimVal(2)-dPixY],'w','LineWidth',3)
    text( xLimVal(2)-lengthScaleBar_pix-dPixX, yLimVal(2)-dPixY-textSpace,[num2str(lengthScaleBar_mu),' \mum'],'HorizontalAlignment','left','color', 'w','FontSize',16)
    % The scale bar for the stresses:
    quiver(xLimVal(2)-lengthScaleBar_pix-dPixX, yLimVal(2)-2*dPixY,fxScaleBar_Pa/forceScale,fyScaleBar_Pa/forceScale,0,'w','LineWidth',2,'MaxHeadSize',5)
    text(  xLimVal(2)-lengthScaleBar_pix-dPixX, yLimVal(2)-2*dPixY-textSpace,[num2str(fxScaleBar_Pa/1000),' kPa'],'HorizontalAlignment','left','color', 'w','FontSize',16)
    hold off
    %!!! Equal axis is important for dimension/y-forces scale to be accuarte!!!
    axis equal
    xlim(xLimVal)
    ylim(yLimVal)
    set(gca,'YDir','reverse','XTick',[],'YTick',[])
    title(['Constrained force field no: ',num2str(frame)])
    
    if doSave==1
        padZeros=3;
        if ~isdir(target_dir)
            mkdir(target_dir)
        end
        filename = [target_dir,filesep,'cellCellForces_',num2str(frame,['%0.',int2str(padZeros),'d'])];
        saveas(gcf,[filename,'.tiff'],'tiffn');
        
        print('-depsc2','-loose', 'frame.eps');    
        str = ['!convert -density 300 -quiet -colorspace rgb -alpha off -depth 8 frame.eps ',filename,'.tif'];
        %str = ['!convert -colorspace rgb -alpha off -depth 8 frame.eps ',filename,'.tif'];
        eval(str);

        % saveas(gcf,[filename, '.eps'], 'psc2');
        % display(['Figure saved to: ',filename,'.tiffn+.eps'])
        %pause(2);        
    end
    pause(1)
    
    if doMovie==1 && (strcmp(movieFormat,'mov') == 1 || strcmp(movieFormat,'MOV'))
        MakeQTMovie('addfigure');
    elseif doMovie==1
        M(k) = getframe;
    end

    if frame~=toDoList(end)
        close(h);
    end
end
!rm frame.eps


if doMovie==1 && (strcmp(movieFormat,'mov') == 1 || strcmp(movieFormat,'MOV'))
    MakeQTMovie('finish');
elseif doMovie==1 && (strcmp(movieFormat,'avi') == 1 || strcmp(movieFormat,'AVI'))
    movie2avi(M,movieFileName,'fps',fps);
end


% Set the position of the current figure to screen size:
% scrsz = get(0,'ScreenSize');
% h     = figure();
% set(h,'Position',scrsz);
% pos   = get(h,'Position');

% figure(h)
% k=1;
% for frame=toDoList
%     axes('Position',[0 0 1 1]);
%     imagesc(I{frame}), title(['Sheet edge of frame: ',num2str(frame)]);
%     hold on
%     plot(sheetBnD(frame).pos(:,2) ,sheetBnD(frame).pos(:,1) ,'k-');
%     plot(sheetEdge(frame).pos(:,2),sheetEdge(frame).pos(:,1),'r.');
%     if frame==toDoList(end)
%         plot(sheetBnD(1).pos(:,2) ,sheetBnD(1).pos(:,1) ,'k-');
%         plot(sheetEdge(1).pos(:,2),sheetEdge(1).pos(:,1),'b.');
%     end
%     colormap gray;
%     set(h,'Position',pos);
%     hold off
%     
% end
