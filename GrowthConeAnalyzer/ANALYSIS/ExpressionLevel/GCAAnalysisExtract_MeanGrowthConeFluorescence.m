function [paramC] = GCAAnalysisExtract_MeanGrowthConeFluorescence(movieData,figuresOn)
%GCAAnalysisExtract_GrowthConeFluorescenceLevel
% INPUT: movieData object
% OUTPUT:
% paramC: Nx1 cell array where N is the mean background subtracted
% intensity value per frame


outDir  = [movieData.outputDirectory_  filesep 'PARAMETER_EXTRACTION' filesep ...
    'Descriptor' filesep 'GrowthCone' filesep 'FluorescenceLevel'];

% if figuresOn == true
%     names{1} = 'BackgroundSub';
%     names{2} = 'GCMask_ForBackEst';
%     names{3} = 'GCMask_ForIntExtract' ;
%     names{4} = 'GCAMask_VeilOnly'; % probably the best estimate of expression correct
%     % how much more
%     %names{5} =
%
%
%
%     for i = 1:numel(names)
%         if ~isdir([movieData.ouptutDirectory_ filesep names{i}]);
%             mkdir([movieData.outputDirectory_ filesep names{i}]);
%         end
%     end
%
% end

% for now skip if file already exists
if exist([outDir filesep 'paramC.mat'],'file')==0;
    
    display(['Extracting Mean Growth Cone Fluorescence For ' movieData.outputDirectory_]);
    if ~isdir(outDir)
        mkdir(outDir)
    end
    load([movieData.outputDirectory_ filesep 'filopodia_fits\Filopodia_Fits_Channel_1' ...
        filesep 'analInfoTestSave.mat']);
    for iFrame = 1:movieData.nFrames_-1
        % load veil est
        img = double(imread([movieData.getChannelPaths{1} filesep movieData.getImageFileNames{1}{iFrame}]));
        %% Make two different growth cone masks
        % One will be slightly more dilated to make sure you don't include too much
        % filo fluorescence in the background estimate - it will use the full
        % length of the entire filo detected via the steerable filter.
        
        maskFiloPreFit = zeros(movieData.imSize_);
        maskFiloAfterFit = zeros(movieData.imSize_);
        % load veil stem mask
        maskVeilStem = analInfo(iFrame).masks.neuriteEdge;
        
        % One will be based on the filo fits and as we will see we might need to
        % dilate more slightly as some of the filo are more than one pixel thick
        
        % load filopodia data for each frame
        filoInfo = analInfo(iFrame).filoInfo;
        
        % Collect all pixels associated with filo
        
        % initiate cell arrays for the filo pixel info
        pixIndicesFiloPreFit = cell(length(filoInfo),1);
        pixIndicesFiloAfterFit = cell(length(filoInfo),1);
        
        % truncate the filopodia by the local fits
        for iFilo = 1:length(filoInfo)
            pixIndicesFiloPreFit{iFilo} = filoInfo(iFilo).('Ext_pixIndices');
            idxEnd = find(pixIndicesFiloPreFit{iFilo} == filoInfo(iFilo).Ext_endpointCoordFitPix);
            pixIndicesFiloAfterFit{iFilo,1} = pixIndicesFiloPreFit{iFilo}(1:idxEnd);
        end
        
        % make conservative mask
        pixIndicesAllAfterFit = vertcat(pixIndicesFiloAfterFit{:});
        % take out NaNs
        pixIndicesAllAfterFit = pixIndicesAllAfterFit(~isnan(pixIndicesAllAfterFit));
        
        % conservative mask for intensity extraction
        maskFiloAfterFit(pixIndicesAllAfterFit) = 1;
        maskFiloAfterFit = logical(maskFiloAfterFit);
        GCMaskStringent = logical(maskFiloAfterFit | maskVeilStem);
        
        %         if figuresOn == true
        %
        %
        %             nx = movieData.imSize_(2);
        %             ny = movieData.imSize_(1);
        %             setFigure(nx,ny,'on');
        %
        %             imshow(img.*GCMaskStringent,[]);
        %             cmap = get(gcf,'colormap');
        %             cmap(1,:) = [0 1 0];
        %             set(gcf,'colormap',cmap);
        %             cDir = [outDir filesep 'TroubleshootFigures' filesep ...
        %                    'GCMaskStringent' ];
        %             if ~isdir(cDir);
        %                 mkdir(cDir);
        %             end
        %                 saveas(gcf,[cDir filesep num2str(iFrame, '%03d') '.png']);
        %                 saveas(gcf,[cDir filesep num2str(iFrame, '%03d') '.fig']);
        %
        %             close gcf
        %         end   % figOn
        
        
        
        
        % make permissive mask
        pixIndicesAllPreFit = vertcat(pixIndicesFiloPreFit{:});
        pixIndicesAllPreFit = pixIndicesAllPreFit(~isnan(pixIndicesAllPreFit));
        maskFiloPreFit(pixIndicesAllPreFit) = 1;
        maskFiloPreFit = logical(maskFiloPreFit);
        GCMaskPermissive = logical(maskFiloPreFit | maskVeilStem);
        
        
        
        
        
        
        
        
        
        % dilate the GCPermissiveMask to get a small mask around the cell
        % for the background subtraction: look at intensity mainly in the
        % direct vicinity surrounding the GCA- specifically did not use a global
        % threshold here as often it misses filopodia and do not want that intensity in the
        % background estiamtion- however do not want to include the whole
        % image as their can be cells etc in the background that you don't
        % want to include as background.
        
        % dilate the mask by 30 pixels to get the area around the growth
        % cone
        GCMaskAreaAroundGC = imdilate(GCMaskPermissive,(strel('disk',30)));
        GCMaskDilSmall = imdilate(GCMaskPermissive,(strel('disk',4))); % just dilate the mask a bit
        GCMaskDilSmall = imfill(GCMaskDilSmall,'holes');
        
        otsuMask = logical(imread([movieData.outputDirectory_ filesep ...
            'masksOtsu' filesep 'maskOtsu' num2str(iFrame,'%03d') '.tif']));
        
        maskBack = ~ GCMaskDilSmall & ~otsuMask ; % get rid of the small dilateed growth cone mask
        % for the background estimate as well as any high intensity blobs
        % that might not be part of the GC but floating in the background
        % estimate area.
        
        maskBack = logical(GCMaskAreaAroundGC.*maskBack); % final background mask just
        % fluorescence information in direct proximity to GC (30 pixels)
        
        % save figure if desired
        if figuresOn == true
            ny = movieData.imSize_(1);
            nx = movieData.imSize_(2);
            setFigure(nx,ny,'off')
            imshow(img.*maskBack,[]);
            cmap = get(gcf,'colormap');
            cmap(1,:) = [0,1,0];
            set(gcf,'colormap',cmap);
            cDir = [outDir filesep 'BackgroundMask'];
            if ~isdir(cDir)
                mkdir(cDir);
            end
            saveas(gcf,[cDir filesep num2str(iFrame,'%03d') '.png']);
            saveas(gcf,[cDir filesep num2str(iFrame, '%03d'),'.eps'],'psc2');
            close(gcf)
        end
        
        backInt(iFrame) =  mean(img(maskBack)); % get the intensity of the background
        
        GCBackSub = img-backInt(iFrame); % subtract the average background intensity for each pixel
        % from the original image to get a background subtacted image
        % clip the mask such that negative values are zero
        GCBackSub(GCBackSub<0) = 0 ;
        paramC{1,iFrame}= mean(GCBackSub(GCMaskStringent));
        
        if figuresOn == 1
            cDir = [outDir filesep 'GCFluorescenceEst'];
            if ~isdir(cDir)
                mkdir(cDir);
            end
            setFigure(nx,ny,'off');
            imshow(GCBackSub.*GCMaskStringent,[]);
            cmap = get(gcf,'colormap');
            cmap(1,:) = [0,1,0];
            set(gcf,'colormap',cmap);
            
            saveas(gcf,[cDir filesep num2str(iFrame,'%03d') '.png']);
            saveas(gcf,[cDir filesep num2str(iFrame, '%03d'),'.eps'],'psc2');
            close gcf
        end
        
    end
    save([outDir filesep 'param_GrowthConeFluorescenceLevel'],'paramC');
else
    display(['Growth Cone Intensities Already Extracted for ' movieData.outputDirectory_ ...
        ' :Skipping'])
end

