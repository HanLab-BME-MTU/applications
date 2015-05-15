function [ output_args ] = GCAVisualsMakeOverlaysFilopodiaMovie(movieData,paramsIn)
%GCAVisualsMakeOverlaysFilopodiaMovie: movieData which overlays the
%segmented filopodia given a certain filopodia filter set
%

plotSpecificFrame =0;
if nargin<2 || isempty(paramsIn)
    paramsIn.colorScheme = 'multi';
end
if nargin<3
    loadFilter = questdlg('Load a filopodia filter set?');
    if strcmpi(loadFilter,'yes')
        filtersRun = searchFiles('filoFilterSet',[], movieData.outputDirectory_, 1,'all',1);
        idxDirectoryInclude  = listSelectGUI(filtersRun,1,'move');
        
        
        paramsIn.filterSetFile = filtersRun{idxDirectoryInclude};
    else
        paramsIn.filterSetFile = [];
    end
    
    %paramsIn.filterSetDirectory = [];  % load filter set
end

if nargin<4
    if ~isempty(paramsIn.filterSetFile)
        % find potential values to plot using this.
        lookInHere =  upDirectory(paramsIn.filterSetFile,1);
        potentialValues = searchFiles('param', [], lookInHere,1,'all',1);
        idx = listSelectGUI(potentialValues,1,'move');
        paramsIn.valuesFile = potentialValues{idx};
        load(paramsIn.valuesFile);
    else
        paramsIn.valuesFile = [];
    end
    
end


nChannels = 1;
for iCh = 1:numel(nChannels)
    imDir = movieData.getChannelPaths{iCh};
    listOfNames = movieData.getImageFileNames{iCh};
    
    %% LOAD SEGMENTATION INFO
    load([movieData.outputDirectory_ filesep 'filopodia_fits' filesep 'Filopodia_Fits_Channel_' num2str(iCh)...
        filesep 'analInfoTestSave.mat']);
    %% LOAD FILOPODIA FILTER INFO (IF DESIRED)
    if ~isempty(paramsIn.filterSetFile)
        % load the filter set
        load(paramsIn.filterSetFile);
    else
        filoFilterSet = []; % no filter.
    end
end
%% Make save directory
if isempty(paramsIn.valuesFile)
    
    saveDir = [movieData.outputDirectory_ filesep 'growth_cone_mask' filesep ...
        'Overlays' paramsIn.colorScheme];
else
    [mainDir,name] = upDirectory(paramsIn.valuesFile,1);
    saveDir = [mainDir filesep name 'Movie'];
end  % if isempty



if ~isdir(saveDir);
    mkdir(saveDir);
end
if plotSpecificFrame == 1
    startFrame = 104;
    endFrame = 104;
else
    startFrame = 1;
    endFrame = numel(listOfNames(:,1))-1;
end

%% Start Plotting Each Frame
for iFrame =startFrame:endFrame
    %% GET THINGS IN ORDER FOR PLOTTING
    
    % get the current filter set
    if ~isempty(filoFilterSet)
        filterSetC = filoFilterSet{iFrame};
    else
        filterSetC = []; %
    end
    
    % Load the veil/stem mask
    veilStemMask = analInfo(iFrame).masks.neuriteEdge;
    roiYXVeilStem = bwboundaries(veilStemMask);
    
    % load the filopodia info for the current frame
    filoInfoC = analInfo(iFrame).filoInfo;
    
    % load the image
    img = double(imread([imDir filesep listOfNames{iFrame}]));
    [ny,nx]= size(img);
    
    if exist('paramC','var');
        valuesC = paramC{iFrame};
    else valuesC = [];
    end
    
    
    %         if rem(ny,2)~=0
    %             % make it even so ffmpeg doesn't flip its shit
    %             ny = ceil(ny+1);
    %         end
    %         if rem(nx,2) ~=0
    %             nx = ceil(nx+1);
    %         end
    
    %% START PLOTTING
    setFigure(nx,ny,'on');
    
    imshow(-img,[])
    
    hold on
    cellfun(@(x) plot(x(:,2),x(:,1),'color','y','Linewidth',2),roiYXVeilStem);
    if strcmpi(paramsIn.colorScheme,'multi');
        
        
        n = length(filoInfoC);
        c = linspecer(n);
        idxRand = randperm(n);
        c = c(idxRand,:);
        count = 1;
        for ifilo = 1:length(filoInfoC)
            
            filterSetIdx = filterSetC(ifilo);
            if filterSetIdx == 1;
                filoInfoIdx = filoInfoC(ifilo);
                valuesIdx = valuesC(count);
                count = count +1;
                
                GCAVisualsFilopodiaOverlaysFromFilterSet(filoInfoIdx,[ny,nx],filterSetIdx,1,c(ifilo,:),valuesIdx,1);
                clear filoInfoIdx filterSetIdx valuesIdx
            end % if filterSetIdx
        end % if ifilo
    else
        GCAVisualsFilopodiaOverlaysFromFilterSet(filoInfoC,[ny,nx],filterSetC,1,[],valuesC); % plot in red
    end
    hold on
    pixelSizeMicron= movieData.pixelSize_/1000;
     forScaleBar = round(10/pixelSizeMicron); % default 10 um scale bar
   
    plotScaleBar(forScaleBar,'Color',[0 0 0],'Location', 'SouthEast');
    
    %  GCAVisualsMakeOverlaysFilopodia(filoInfoC,[ny,nx],1,1,[],0);
    saveas(gcf,[saveDir filesep num2str(iFrame,'%03d') '.png']);
    %saveas(gcf,[saveDir filesep 'FilopodiaOverlay' num2str(iFrame,3) '.eps'],'psc2');
    saveas(gcf,[saveDir filesep num2str(iFrame,'%03d')]);
    close gcf
    
end % iFrame
%
cd(['C:\Users\Maria\Desktop\' ...
    'ffmpeg-20141210-git-ae81680-win64-static\ffmpeg-20141210-git-ae81680-win64-static\bin']);
[~,num] = upDirectory(movieData.outputDirectory_,2,1);
[~,date] = upDirectory(movieData.outputDirectory_,3,1);
[~,cond ]= upDirectory(movieData.outputDirectory_,4,1);
execute = ['ffmpeg -y -r 5 -i ' saveDir filesep '%03d.png' ...
    ' -crf 22 -pix_fmt yuv420p -b 20000k -bt 20000k ' saveDir filesep [cond date num] '.mp4'];
system(execute)
cd(saveDir)

end

