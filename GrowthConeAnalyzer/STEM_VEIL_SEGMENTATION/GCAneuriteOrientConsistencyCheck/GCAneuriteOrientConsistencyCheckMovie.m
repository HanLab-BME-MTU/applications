function [ backboneInfo,frames2Fix] = GCAneuriteOrientConsistencyCheckMovie(movieData,paramsIn)
%% GCAneuriteOrientConsistencyCheckMovie.m (Was checkForConsistencyNeuriteOrientMovie until 20140529)
% Note need to fix the input such that this is more formally an appropriate
% movieWrapperFunction

% This function inputs the backboneInfo.mat data structure output from
% getNeuriteOrientationMovie.m and tests for consistency among the backbone
% seed masks. It makes two main assumptions
% (1) that the neurite entrance point orientation is relatively
% constant throughout the whole movie
% (2) that the majority of backbone seed masks were of the correct
% orientation, so that any outliers found need to be corrected-
% the next best neurite seed candidate that complies with the majority
% vote is chosen
%


% Input:
%
%   movieData - The MovieData object describing the movie, as created using
%   setupMovieDataGUI.m
%
%   paramsIn - Structure with inputs for optional parameters. The
%   parameters should be stored as fields in the structure, with the field
%   names and possible values as described below:
%
%   Possible Parameter Structure Field Names:
%       ('FieldName' -> possible values)
%
%       ('OutputDirectory' -> character string) Optional. A character
%       string specifying the directory to save the backboneInfo structure to.
%       If not input, the backboneInfo will be saved in the same directory
%       as the movieData, in a sub-directory called "neurite_orientation_estimation_fixes"
%
%       ('ChannelIndex'-> Positive integer scalar or vector) The integer
%       indices of the channel(s) on which to perform the neurite orientation est.
%       This index corresponds to the channel's location in the array
%       movieData.channels_. If not input, all channels will be analyzed
%
%       ('ProcessIndex' -> Positive integer scalar or vector)
%       This parameter specifies the output process to use for performing the
%       estimation
%       This allows a previously processed image (ie shade corrected or
%       background subtracted to potentially be input). If not input, the
%       backbone information will be calculated from the channels (ie raw
%       images)
%
%       ('SizeOfConsistencyRestraint' -> scalar)
%       Default = 10 pixels
%       This parameter allows the user to control the amount of noise is
%       allowed for the consistency metric. The function currently identifies
%       outliers by simply dilating a small disk mask around the backbone
%       seed entrance points. backbone seed entrance coordinates that fall
%       outside of the resulting connected component are considered
%       potential outliers to be corrected. Smaller values for this parameter
%       are used when the true neurite entrance varies only slightly throughout
%       the movie, while larger values are sometimes required for neurites
%       backbones that are very mobile.
%
%       ('Troubleshoot Plots' -> if true will make three troubleshoot
%       folders with the larger outputDirectory for the given channel
%

%%

if nargin < 1 || ~isa(movieData,'MovieData')
    error('The first input must be a valid MovieData object!')
end

if nargin < 2
    paramsIn.OutputDirectory = [movieData.outputDirectory_ filesep 'neurite_orientation_estimation_fixes'];
    paramsIn.ChannelIndex = 1;
    paramsIn.ProcessIndex = 0; % use raw images
    paramsIn.SizeOfConsistencyRestraint = 10;% default 10 
    paramsIn.plots = 1;
end


% FOR WHEN MAKE PROCESS
%Get the indices of any previous mask refinement processes from this function
% iProc = movieData.getProcessIndex('GetNeuriteOrientationProcess',1,0);
%
% %If the process doesn't exist, create it
% if isempty(iProc)
%     iProc = numel(movieData.processes_)+1;
%     movieData.addProcess(GetNeuriteOrientationProcess(movieData,movieData.outputDirectory_));
% end
%
% %Parse input, store in parameter structure
% p = parseProcessParams(movieData.processes_{iProc},paramsIn);

p = paramsIn;


%% Init:
nFrames = movieData.nFrames_;
nChan = numel(p.ChannelIndex);
imSize = movieData.imSize_;
ySize = imSize(1);
xSize = imSize(2);
mkClrDir(p.OutputDirectory)




%% Start Fixing Backbone Orientations

for iCh = 1:nChan
    display(['Testing For Neurite Orientation Consistency Channel ' num2str(iCh)]);
    
    % load backboneInfo  %%% NEED TO CHANAGE WHEN MAKE A PROCESS
    load([movieData.outputDirectory_ filesep 'neurite_orientation_estimations' filesep ...
        'Neurite_Backbone_Seed_Channel_' num2str(iCh) filesep 'backboneInfo.mat'] );
    
    % make final output dir where backboneInfo will be saved
    saveDir =  [p.OutputDirectory filesep 'Neurite_Backbone_Seed_Channel_' num2str(iCh)];
    saveDirOld = [movieData.outputDirectory_ filesep 'neurite_orientation_estimations' filesep 'Neurite_Backbone_Seed_Channel_' num2str(p.ChannelIndex)]; % transfer over some files from
    % old files
    
    mkClrDir(saveDir)
    
    % get the list of image filenames
    if p.ProcessIndex == 0
        imDir = movieData.channels_(iCh).channelPath_;
    else
        imDir = movieData.proceses_{p.ProcessIndex}.outfilePaths_;
    end
    
    listOfImages = searchFiles('.tif',[],imDir,0);
    if isempty(listOfImages)
        error('No Images Found: Check Input Directory');
    end
    
    % if restart == 1
    %    % startFrame = numel(analInfo)-1;
    %
    % else
    %     startFrame = 1;
    % end
    
    %
    
    
    
    % get all the coordsin for all the backbone seeds for all movies
    
    coordIn = arrayfun(@(x) x.coordsEnterNeurite,backboneInfo,'uniformoutput',0);
    
    coordInPix = cellfun(@(i) sub2ind([ySize,xSize],i(2),i(1)), coordIn,'uniformoutput',0);
    % put into a mask
    dilateMask = zeros([ySize,xSize]);
    dilateMask(vertcat(coordInPix{:}))=1;
    
    dilateMask = imdilate(dilateMask,strel('disk',10)); % changed to 10
    
    
    CCDil  = bwconncomp(dilateMask);
    frames2Fix = [];
    if(CCDil.NumObjects>1)
        allPixels = vertcat(coordInPix{:});
        idxPixelsInClust = cellfun(@(x) intersect(allPixels,x),CCDil.PixelIdxList,'uniformoutput',0);
        for iCluster = 1:numel(idxPixelsInClust)
            current = idxPixelsInClust{iCluster};
            numTimesInClust(iCluster) =  sum(arrayfun(@(i) sum(allPixels==i),current)) ;
        end
        % count per frame
        % numPixInClust = cellfun(@(x) length(intersect(allPixels,x)),CCDil.PixelIdxList);
        % pixelsClustLarge = CCDil.PixelIdxList{numPixInClust==max(numPixInClust)};
        pixelsClustLarge = vertcat(idxPixelsInClust{numTimesInClust==max(numTimesInClust)});
        pixelsNOClust = vertcat(idxPixelsInClust{numTimesInClust~=max(numTimesInClust)});
        
        test = mat2cell(pixelsNOClust,ones(length(pixelsNOClust),1));
        frames2Fix =  cellfun(@(x) find(allPixels==x),test,'uniformoutput',0);
        frames2Fix = vertcat(frames2Fix{:});
        
    end
    frames2Fix = sort(frames2Fix);
    if isempty(frames2Fix)
        display('All Neurite Orientations Consistent: No Changes Were Made')
        
    else
        display('Fixing Minority Neurite Orientations to the Majority')
        %% Troubleshooting plots
        if p.plots == 1
            % copy all bad files into a new dir
            badFrameFile = [saveDir filesep 'Frames2FixBefore'];
            if ~isdir(badFrameFile)
                mkdir(badFrameFile)
            end
            
            name{1} = 'BeforeAndAfterConnect';
            name{2} = 'CandSeeds';
            name{3} = 'RidgeCandBeforeAfterClean';
            
            
            for iFrame = 1:length(frames2Fix)
                
                for iTransfer = 1:3
                    if isdir([saveDirOld filesep name{iTransfer}])
                        source = [saveDirOld filesep name{iTransfer} filesep num2str(frames2Fix(iFrame),'%03d') '.tif'];
                        
                        dest = [badFrameFile filesep name{iTransfer} num2str(frames2Fix(iFrame),'%03d') '.tif'] ;
                        copyfile(source,dest);
                        
                        
                    else
                        warning( [ name{iTransfer} 'troubleshoot plots are copied from the previous step- no old troubleshoot plots found'])
                    end % isdir
                end % iTransfer
            end % iFrame
        end % if p.plots
        %% fix backbone seed by finding the ridge closest to the majority cluster
        % make a mask of just the majority cluster
        maskMClust = zeros(ySize,xSize);
        maskMClust(pixelsClustLarge) = 1;
        % dist transform of that mask
        distTransFromMClust = bwdist(maskMClust);
        
        
        
        %%
        for iFrame = 1:length(frames2Fix)
            display(['Fixing Frame ' num2str(frames2Fix(iFrame))])
            % load the candidate ridgeMask (mask after cleaning and linear connectivity)
            candRidges =  backboneInfo(frames2Fix(iFrame)).bodyReconstruct.AfterConnect;
            
            %
            origBBMask = backboneInfo(frames2Fix(iFrame)).backboneSeedMask;
            % set the original backbone to 0
            candRidges(origBBMask==1) = 0;
            
            
            
            % get the connected component structure of the mask
            CCNMS = bwconncomp(candRidges);
            % get their labels
            labelsRidges = labelmatrix(CCNMS);
            
            % apply the ridge mask
            distRidgeMat = candRidges.*distTransFromMClust;
            %    distances = distances(:);
            %    distances = distances(distances~=0);
            %    minDist = min(distances);
            % find all ridges within 50 pixels
            candLabels = labelsRidges(distRidgeMat<20&distRidgeMat~=0);
            candLabels = unique(candLabels);
            % get the majority body
            % collect backbone from all good frame
            goodFrames = setdiff(1:nFrames,frames2Fix);
            backbones = arrayfun(@(x) x.backboneSeedMask,backboneInfo(goodFrames),'uniformoutput',0);
            sumBB = zeros(ySize,xSize);
            for iBB = 1:length(goodFrames)
                sumBB = sumBB+backbones{iBB}; % likely better way than a loop for this but quick fix
            end
            % find pixels in structure that were more in more than 5 frames (to remove
            % outliers)
            sumBB(sumBB<=5) =0;
            sumBB(sumBB>0) = 1;
            pixMajBB = find(sumBB==1) ;
            % find max overlap between cand ridge and
            overlap = cellfun(@(x) length(intersect(pixMajBB,x)),CCNMS.PixelIdxList(candLabels));
            labelKeep = candLabels(overlap==max(overlap));
            %% try to make flag  so if the wrong side was connected can fix. (as in DOCK01 new data)
            if isempty(labelKeep) % sometimes it connects to the wrong side
                % test the original candidate to see if there is significant overlap
                % sometimes simply have a long candidate attach to wrong side
                overlapOrig = length(intersect(pixMajBB,find(origBBMask==1))) ;
                
                % flag to fix
                
                
                
                
                
                if ~isempty(overlapOrig)
                    % get other endpoint
                    EPCoordsBodyMask = getEndpoints(find(origBBMask==1),[ySize,xSize]);
                    idxEPs = sub2ind([ySize,xSize],EPCoordsBodyMask(:,2),EPCoordsBodyMask(:,1));
                    enterCoords = backboneInfo(frames2Fix(iFrame)).coordsEnterNeurite;
                    enterIdx  = sub2ind([ySize,xSize],enterCoords(:,2),enterCoords(:,1));
                    idxEPs(idxEPs ==enterIdx) = [];
                    pixBackboneNew = find(origBBMask==1) ;
                    useOtherSide =1 ;
                end
                %      end
            else
                pixBackboneNew = CCNMS.PixelIdxList{labelKeep};
                useOtherSide =0;
            end
            
            
            %    candRidgeIdx = find(cellfun(@(x) ~isempty(intersect(x,pixelsClustLarge)),CCNMS.PixelIdxList));
            %    toTest = CCNMS.PixelIdxList(candRidgeIdx) ;
            %    sizeRidge = cellfun(@(x) length(x),toTest);
            %    candRidgeIdx = candRidgeIdx(sizeRidge==max(sizeRidge));
            %
            %
            
            
            backboneSeed = zeros(ySize,xSize);
            
            
            
            
            backboneSeed(pixBackboneNew)=1;
            
            
            
            
            
            dims = [ySize,xSize];
            boundaryMask(1:dims(1),1) =1;
            boundaryMask(1:dims(1),dims(2))=1;
            boundaryMask(1,1:dims(2))= 1;
            boundaryMask(dims(1),1:dims(2)) =1;
            idxEnterNeurite = find(backboneSeed ==1 & boundaryMask ==1);
            [yEnter,xEnter] = ind2sub([ySize,xSize],idxEnterNeurite);
            if isempty(idxEnterNeurite) % need to interpolate to make sure closed contour
                % interpolate between to nearest point on boundary
                % find endpoint of backboneSeed
                
                EPsBackbone = getEndpoints(pixBackboneNew,[ySize,xSize]);
                if useOtherSide == 1; % flag to not use old endpoint
                    idxEPs = sub2ind([ySize,xSize],EPsBackbone(:,2),EPsBackbone(:,1));
                    enterCoords = backboneInfo(frames2Fix(iFrame)).bodyReconstruct.coordsEnterNeurite;
                    enterIdx  = sub2ind([ySize,xSize],enterCoords(:,2),enterCoords(:,1));
                    idxEPs(idxEPs ==enterIdx) = [];
                    [EPY,EPX] = ind2sub([ySize,xSize],idxEPs);
                    EPsBackbone = [EPX, EPY];
                end
                
                [yBoundary,xBoundary] = find(boundaryMask==1);
                [idxBoundaryClose,dist] = KDTreeBallQuery([xBoundary,yBoundary] ,EPsBackbone,20);
                % quick fix
                if isempty(idxBoundaryClose{1})
                    [idxBoundaryClose,dist] = KDTreeBallQuery([xBoundary,yBoundary],EPsBackbone,100);
                end
                distAll =  vertcat(dist{:});
                
                % find the minimum distance
                toSave = find(distAll ==min(distAll));
                if length(toSave) > 1
                    toSave = toSave(1);
                end
                E = arrayfun(@(i) [repmat(i, [numel(idxBoundaryClose{i}) 1]) idxBoundaryClose{i}], 1:length(EPsBackbone(:,1)), 'UniformOutput', false);
                E = vertcat(E{:});
                linkCoords = bresenham([EPsBackbone(E(toSave,1),1),EPsBackbone(E(toSave,1),2)],[xBoundary(E(toSave,2)),yBoundary(E(toSave,2))]);
                backboneSeed(linkCoords(:,2), linkCoords(:,1) ) = 1;
                linked = 1;
                xEnter = xBoundary(E(toSave,2));
                yEnter = yBoundary(E(toSave,2));
                clear linkedCoords E
            end    % need to close contour
            % save new neurite entrance and the new backboneSeed
            backboneInfo(frames2Fix(iFrame)).backboneSeedMask= backboneSeed;
            backboneInfo(frames2Fix(iFrame)).coordsEnterNeurite = [xEnter,yEnter];
            backboneInfo(frames2Fix(iFrame)).timeStamp = clock;
            % document in a folder
            close gcf
            %% Troubleshoot Plots
            if p.plots  == 1
                fixDir = [saveDir filesep 'AfterFix' ];
                
                if ~isdir(fixDir)
                    mkdir(fixDir)
                end
                
                img = double(imread([listOfImages{frames2Fix(iFrame),2} filesep listOfImages{frames2Fix(iFrame),1}]));
                setFigure(xSize,ySize,'off')
                
                imshow(img,[]) ;
                hold on
                spy(origBBMask,'g')
                hold on
                spy(backboneSeed,'r');
                saveas(gcf,[fixDir filesep 'OldVsNew' num2str(frames2Fix(iFrame),'%03d') '.tif']);
            end
            close gcf
            
        end
        
        
    end
      
        save([saveDir filesep 'backboneInfoFix.mat'],'backboneInfo');
        save([saveDir filesep 'framesFixed.mat'],'frames2Fix');
        save([saveDir filesep 'paramsIn.mat'],'paramsIn');
end % iCh
end % The END

