function checkCutOutCells(constrForceField,forceField,imageFileList,frame)
load('fileAndFolderNames.mat')

if nargin < 1 || isempty(constrForceField)
    try
        load(path_cellCellForces)
    catch
        [filename, pathname] = uigetfile({'*.mat';'*.*'}, ...
            'Select cellCellForces.mat to be used');
        fileStruct=load([pathname filesep filename]);
        constrForceField=fileStruct.constrForceField;
    end
end

if nargin < 2 || isempty(forceField)
    try
        load(path_forceField)
    catch
        [filename, pathname] = uigetfile({'*.mat';'*.*'}, ...
            'Select forceField.mat to be used');
        fileStruct=load([pathname filesep filename]);
        forceField=fileStruct.forceField;
    end
end

%read in stack of image files:
if nargin < 3 || isempty(imageFileList)
    try
        imageFileList=getFileListFromFolder(path_XtraFinal);
    catch
        [filename, pathname] = uigetfile({'*.tif';'*.TIF';'*.*'}, ...
            'Select the first image file to be overlayed');
        if ~ischar(filename) || ~ischar(pathname)
            return;
        end
        imageFileList = getFileStackNames([pathname filesep filename]);
    end
else
    isValid = 1;
    for iframe = 1:numel(imageFileList)
        isValid = isValid && exist(imageFileList{iframe}, 'file');
    end
    if ~isValid
        error('Invalid input files.');
    end
end


toDoList=[];
if nargin < 4 || isempty(frame)
    toDoList=input('Which frame should be analyzed [all]: ');
    if isempty(toDoList)
        for frame=1:length(constrForceField)
            if ~isempty(constrForceField{frame})
                toDoList=horzcat(toDoList,frame);
            end
        end
    end
else
    toDoList=frame;
end
    
display(['Do frames: ',num2str(toDoList)]);

% first test if the cluster analysis has been performed on this frame. For
% example, there is only a single cell:
% if ~isfield(constrForceField{frame},'clusterAnalysis')
%     display(['For frame: ',num2str(frame),' no cluster analysis has been performed!'])
%     display('Nothing to do!')
%     return;
% end
for frame =toDoList
    currentImage = double(imread(imageFileList{frame}));
    dPix=50;
    min_x=min(constrForceField{frame}.roi.pos(:,1))-dPix;
    max_x=max(constrForceField{frame}.roi.pos(:,1))+dPix;
    min_y=min(constrForceField{frame}.roi.pos(:,2))-dPix;
    max_y=max(constrForceField{frame}.roi.pos(:,2))+dPix;
    maxForcePlot=1/forceField(frame).par.gridSpacing*max(sqrt(forceField(frame).vec(:,1).^2+forceField(frame).vec(:,2).^2));
    [rows,cols]=size(currentImage);
    
    factor_Pa2nN=constrForceField{frame}.par.gridSpacing^2*constrForceField{frame}.par.pixSize_mu^2*10^(-3);
    subplot(1,2,1)
    imagesc(currentImage)
    colormap('gray')
    hold on
    %plot innter boundary:
    plot(constrForceField{frame}.segmRes.curve(:,1),constrForceField{frame}.segmRes.curve(:,2),'k')
    %plot complete forceField:
    quiver(forceField(frame).pos(:,1),forceField(frame).pos(:,2),forceField(frame).vec(:,1)/maxForcePlot,forceField(frame).vec(:,2)/maxForcePlot,0,'g');
    marker=['r','b','m','c','g','y','k'];
    if isfield(constrForceField{frame},'cell')
        for k=1:length(constrForceField{frame}.cell)
            %plot cell{j}:
            quiver(constrForceField{frame}.cell{k}.pos(:,1),constrForceField{frame}.cell{k}.pos(:,2),constrForceField{frame}.cell{k}.vec(:,1)/maxForcePlot,constrForceField{frame}.cell{k}.vec(:,2)/maxForcePlot,0,marker(mod(k,7)+1));
            %     if j==1 % if there are only two cells plot the residual forces at the interface:
            %         resForcePos=constrForceField{frame}.cell{k}.interface.center;
            %     else % if there are more than two cells plot the residual forces at center of mass of each cell:
            %         resForcePos=constrForceField{frame}.cell{k}.stats.resForce.pos;
            %     end
            %!!!            % convert the force back to a stress to plot it with the
            % traction stress. That is a bit dirty!
            % quiver(resForcePos(1),resForcePos(2),constrForceField{frame}.cell{k}.stats.resForce.vec(:,1)/(maxForcePlot*factor_Pa2nN),constrForceField{frame}.cell{k}.stats.resForce.vec(:,2)/(maxForcePlot*factor_Pa2nN),0,marker(mod(k,7)+1));
            % plot(resForcePos(1),resForcePos(2),['o',marker(mod(k,7)+1)])
            %plot boundaries:
            plot(constrForceField{frame}.cell{k}.boundary(:,1),constrForceField{frame}.cell{k}.boundary(:,2),marker(mod(k,7)+1))
        end
        for k=1:length(constrForceField{frame}.interface)
            plot(constrForceField{frame}.interface{k}.pos(:,1),constrForceField{frame}.interface{k}.pos(:,2),'--w')
        end
    end
    % plot the holes:
    if isfield(constrForceField{frame}.segmRes,'hole')
        for holeId=1:length(constrForceField{frame}.segmRes.hole)
            plot(constrForceField{frame}.segmRes.hole{holeId}.curve(:,1),constrForceField{frame}.segmRes.hole{holeId}.curve(:,2),'k')
        end
    end
    % plot the error for this cluster. Use the same conversion
    % factor factor_Pa2nN as for interfacial stresses:
    quiver(min_x+2*dPix,min_y+2*dPix,constrForceField{frame}.errorSumForce.vec(1)/(maxForcePlot*factor_Pa2nN),constrForceField{frame}.errorSumForce.vec(2)/(maxForcePlot*factor_Pa2nN),0,'w');
    hold off
    title(['Constrained force field no: ',num2str(frame)])
    axis equal
    xlim([max([1 min_x]) min([cols,max_x])])
    ylim([max([1 min_y]) min([rows,max_y])])
    set(gca,'YDir','reverse')
    
    subplot(1,2,2)
    imagesc(currentImage)
    title(['Plain Ecad image, frame: ',num2str(frame)])
    axis equal
    xlim([max([1 min_x]) min([cols,max_x])])
    ylim([max([1 min_y]) min([rows,max_y])])
    set(gca,'YDir','reverse')
    
    
    input(['Hit enter to continue [Frame= ',num2str(frame),';  Rd= ',num2str(constrForceField{frame}.segmRes.params.dilationR),';  errF= ',num2str(round(constrForceField{frame}.errorSumForce.mag)),'nN]:'])
end