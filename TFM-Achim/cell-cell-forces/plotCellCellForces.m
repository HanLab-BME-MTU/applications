function plotCellCellForces(constrForceField,forceField,toDoList,imageFileList,target_dir,doSave)
% This function plots the residual forces in the cell cluster!


if nargin < 1 || isempty(constrForceField)
   [filename_constrForceField, pathname_constrForceField] = uigetfile({'*.mat';'*.*'}, ...
       'Select cellCellForces.mat to be used');
       %the stage drift Transformation:
       fileStruct=load([pathname_constrForceField filesep filename_constrForceField]);
       constrForceField=fileStruct.constrForceField;
       
       target_dir=pathname_constrForceField;
end

if nargin < 2 || isempty(forceField) 
   [filename_forceField, pathname_forceField] = uigetfile({'*.mat';'*.*'}, ...
       'Select forceField.mat to be used');
       %the stage drift Transformation:
       fileStruct=load([pathname_forceField filesep filename_forceField]);
       forceField=fileStruct.forceField;
end

if nargin<3 || isempty(toDoList)
    toDoList=[];
    for frame=1:length(constrForceField)
        if ~isempty(constrForceField{frame}) && isfield(constrForceField{frame},'cell')
            toDoList=horzcat(toDoList,frame);
        end
    end
end

%read in Stack of images:
if nargin < 4 || isempty(imageFileList)
   [filename, pathname] = uigetfile({'*.tif';'*.jpg';'*.png';'*.*'}, ...
       'Select First Ecad-Image');
   
   if ~ischar(filename) || ~ischar(pathname)
       return;
   end
   
   imageFileList = getFileStackNames([pathname filesep filename]);
elseif isdir(imageFileList)
    imageFileList=getFileListFromFolder(imageFileList);
else
    isValid = 1;
    for i = 1:numel(imageFileList)
        isValid = isValid && exist(imageFileList{frame}, 'file');
    end
    if ~isValid
        error('Invalid input files.');
    end
end

%get the target directory:
if nargin < 5 || isempty(target_dir)
    target_dir = uigetdir('','Please select target directory');
end

if nargin < 6 || isempty(doSave)
    doSave = 0;
end


%% Determine the plot boundaries and pos of scale bars.
dPix=125;
dPixX=100;
dPixY=50;
textSpace=30;
    
currentImage = double(imread(imageFileList{toDoList(1)}));
[rows,cols]=size(currentImage);
xmin=Inf;
xmax=0;
ymin=Inf;
ymax=0;
forceScale=0;
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
end
xLimVal=[max([1 xmin]) min([cols,xmax])];
yLimVal=[max([1 ymin]) min([rows,ymax])];

% set the values for the scale bars:
lengthScaleBar_mu=5;
fxScaleBar_Pa=15000;
fyScaleBar_Pa=0;


%% plot:
for frame=toDoList
    currentImage = double(imread(imageFileList{frame}));
    
    % calculate the size of um in pixel for a scale bar:
    pixSize_mu=forceField(frame).par.pixSize_mu;
    lengthScaleBar_pix=lengthScaleBar_mu/pixSize_mu;    

    
    % Number of cells and interfaces to plot:
    numCells=length(constrForceField{frame}.cell);
    numIntf =length(constrForceField{frame}.interface);
    
    % scaling factor to convert residual forces [nN] back to stresses:
    factor_Pa2nN=constrForceField{frame}.par.gridSpacing^2*constrForceField{frame}.par.pixSize_mu^2*10^(-3);
    
    % Here starts the image:
    scrsz = get(0,'ScreenSize');
    h=figure('Position',[1 scrsz(4) scrsz(3) scrsz(4)]);
    %figure(frame)
    imagesc(currentImage)
    colormap('gray')
    hold on
    %plot innter boundary:
    plot(constrForceField{frame}.segmRes.curve(:,1),constrForceField{frame}.segmRes.curve(:,2),'k')
    %plot complete forceField:
    quiver(forceField(frame).pos(:,1),forceField(frame).pos(:,2),forceField(frame).vec(:,1)/forceScale,forceField(frame).vec(:,2)/forceScale,0,'g');
    marker=['r','b','m','c','g','y','k'];
    for k=1:numCells
        %plot cell{j}:
        quiver(constrForceField{frame}.cell{k}.pos(:,1),constrForceField{frame}.cell{k}.pos(:,2),constrForceField{frame}.cell{k}.vec(:,1)/forceScale,constrForceField{frame}.cell{k}.vec(:,2)/forceScale,0,marker(mod(k,7)+1));
%!!!    % convert the force back to a stress to plot it with the
        % traction stress. That is a bit dirty!
        resForcePos=constrForceField{frame}.cell{k}.stats.resForce.pos;
        quiver(resForcePos(1),resForcePos(2),constrForceField{frame}.cell{k}.stats.resForce.vec(:,1)/(forceScale*factor_Pa2nN),constrForceField{frame}.cell{k}.stats.resForce.vec(:,2)/(forceScale*factor_Pa2nN),0,marker(mod(k,7)+1));
        plot(resForcePos(1),resForcePos(2),['o',marker(mod(k,7)+1)])
        %plot boundaries:
        plot(constrForceField{frame}.cell{k}.boundary(:,1),constrForceField{frame}.cell{k}.boundary(:,2),marker(mod(k,7)+1))
    end
    for k=1:numIntf % there are only 'j' interfaces:
        plot(constrForceField{frame}.interface{k}.pos(:,1),constrForceField{frame}.interface{k}.pos(:,2),'--w')
    end
    % plot the holes:
    if isfield(constrForceField{frame}.segmRes,'hole')
        for holeId=1:length(constrForceField{frame}.segmRes.hole)
            plot(constrForceField{frame}.segmRes.hole{holeId}.curve(:,1),constrForceField{frame}.segmRes.hole{holeId}.curve(:,2),'k')
        end
    end
    % The scale bar um/pix:
    plot([xLimVal(2)-lengthScaleBar_pix-dPixX  xLimVal(2)-dPixX], [yLimVal(2)-dPixY yLimVal(2)-dPixY],'w','LineWidth',3)
    text( xLimVal(2)-lengthScaleBar_pix-dPixX, yLimVal(2)-dPixY-textSpace,[num2str(lengthScaleBar_mu),' \mum'],'HorizontalAlignment','left','color', 'w','FontSize',16)
    % The scale bar for the stresses:
    quiver(xLimVal(2)-lengthScaleBar_pix-dPixX, yLimVal(2)-2*dPixY,fxScaleBar_Pa/forceScale,fyScaleBar_Pa/forceScale,0,'w','LineWidth',2,'MaxHeadSize',5)
    text(  xLimVal(2)-lengthScaleBar_pix-dPixX, yLimVal(2)-2*dPixY-textSpace,[num2str(fxScaleBar_Pa/1000),' kPa'],'HorizontalAlignment','left','color', 'w','FontSize',16)
    hold off
    axis equal
    xlim(xLimVal)
    ylim(yLimVal)
    set(gca,'YDir','reverse')
    title(['Constrained force field no: ',num2str(frame)])

    if doSave==1
        padZeros=3;
        if ~isdir(target_dir)
            mkdir(target_dir)
        end
        filename = [target_dir,filesep,'cellCellForces_',num2str(frame,['%0.',int2str(padZeros),'d'])];
        saveas(gcf,[filename,'.tiff'],'tiffn');
        %saveas(gcf,[filename, '.eps'], 'psc2');
        display(['Figure saved to: ',filename,'.tiffn+.eps'])
        %pause(2);        
    end
    pause(1)
    if frame~=toDoList(end)
        close(h);
    end
end