function [ backboneInfoFix,frames2Fix] = GCAneuriteOrientConsistencyCheckMovie(movieData,varargin)
%% GCAneuriteOrientConsistencyCheckMovie: MovieDataWrapper for
%  GCAneuriteOrientConsistencyCheck: (STEP II of GCA Segmentation)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUT:
%
%   movieData - The MovieData object describing the movie, as created using
%   setupMovieDataGUI.m
%
%   Possible Parameter Structure Field Names:
%       ('FieldName' -> possible values)
%
% Generic Fields:
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
% 
%
% Specific Input:
%       See GCAgetNeuriteOrientConsistencyCheck.m for details.
%
%% INPUT PARSER
if nargin < 1 || ~isa(movieData,'MovieData')
    error('The first input must be a valid MovieData object!')
end

%%Input check
ip = inputParser;

ip.CaseSensitive = false;

% PARAMETERS
% Generic
defaultOutDir = [movieData.outputDirectory_ filesep...
    'SegmentationPackage' filesep 'StepsToReconstruct' filesep 'II_neurite_orientation_refinements'];

ip.addParameter('OutputDirectory',defaultOutDir,@(x) ischar(x));
ip.addParameter('ChannelIndex',1);
ip.addParameter('ProcessIndex',0);

% Specific
ip.addParameter('TSOverlays',true,@(x) islogical(x));
ip.addParameter('SizeOfConsistencyRestraint',5,@(x) isscalar(x));
ip.addParameter('CheckOrient',false,@(x) islogical(x));

ip.parse(varargin{:});
p = ip.Results;
%%
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
%%
%% Init:
channels = ip.Results.ChannelIndex;
nChan = numel(ip.Results.ChannelIndex);
sizeImg = movieData.imSize_;
xSize = sizeImg(2); 
ySize= sizeImg(1); 
%%
for iCh = 1:nChan
    
    display(['Testing For Neurite Orientation Consistency Channel ' num2str(channels(iCh))]);
    
    saveDir = [ip.Results.OutputDirectory filesep 'Channel_' num2str(channels(iCh))];
    
    if ~isdir(saveDir)
        mkdir(saveDir)
    end
    
    % load backboneInfo  %%% NEED TO CHANAGE WHEN MAKE A PROCESS: check
    % previous process - load and make sure run through completely
    
    saveDirOld = [movieData.outputDirectory_ filesep ...
        'SegmentationPackage' filesep 'StepsToReconstruct' filesep ...
        'I_neurite_orientation' filesep 'Channel_' num2str(iCh) ];
    load([saveDirOld filesep 'backboneInfo.mat']);
    
    
    
    % old files
    
    if isdir(saveDir)
        mkdir(saveDir)
    end
    
    % get the list of image filenames
    if p.ProcessIndex == 0
        imDir = movieData.channels_(iCh).channelPath_;
    else
        imDir = movieData.proceses_{ip.Results.ProcessIndex}.outfilePaths_;
    end
    
    listOfImages = searchFiles('.tif',[],imDir,0);
    if isempty(listOfImages)
        error('No Images Found: Check Input Directory');
    end
    %% Run Function
    if ip.Results.CheckOrient == true; % initiate the figure
       setFigure(xSize,ySize,'on'); 
%     
    img = double(imread([movieData.getChannelPaths{1} filesep movieData.getImageFileNames{1}{1}]));
    imshow(-img,[]);
    hold on
    end  
    
    [backboneInfoFix,frames2Fix,modified]=  GCAneuriteOrientConsistencyCheck(backboneInfo,p);
   
    p.modified = modified;
    %% Troubleshoot Plots
    if (ip.Results.TSOverlays  == true && ~ isempty(frames2Fix))
        fixDir = [saveDir filesep 'AfterFix' ];
        
        if ~isdir(fixDir)
            mkdir(fixDir)
        end
        
        for iFrame = 1:length(frames2Fix)
            
            img = double(imread([listOfImages{frames2Fix(iFrame),2} filesep listOfImages{frames2Fix(iFrame),1}]));
            setFigure(xSize,ySize,'on');
            
            imshow(-img,[]) ;
            hold on
            % plot the original signal 
            spy(backboneInfo(frames2Fix(iFrame)).linkedRidgesFinal,'b');  
            
            % plot the alignment mask 
            spy(backboneInfoFix(frames2Fix(iFrame)).alignmentMask,'m');
            
            % plot the origBB 
            origBBMask = backboneInfo(frames2Fix(iFrame)).backboneSeedMask;
            spy(origBBMask,'g')
            hold on
            
            % plot the final BB seed
            backboneSeed = backboneInfoFix(frames2Fix(iFrame)).backboneSeedMask;
            spy(backboneSeed,'r');
            
            % plot the old and new input neurite coords 
            [coordsOrg] = backboneInfo(frames2Fix(iFrame)).coordsEnterNeurite;
            scatter(coordsOrg(1),coordsOrg(2),'g','filled');
            [coordsNew] = backboneInfoFix(frames2Fix(iFrame)).coordsEnterNeurite ;
            scatter(coordsNew(1),coordsNew(2),'y','filled');
           
            text(10,10,{'Yellow Marks' ; 'Corrected Entry Point'}, 'Color','k');
            
            saveas(gcf,[fixDir filesep 'OldVsNew' num2str(frames2Fix(iFrame),'%03d') '.tif']);
            close gcf
        end % iFrame
        
        % copy some of the original files over.
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
    end % if ip.Results
    %% Save Information
    backboneInfo  = backboneInfoFix; 
    save([saveDir filesep 'backboneInfoFix.mat'],'backboneInfo');
    save([saveDir filesep 'framesFixed.mat'],'frames2Fix');
    save([saveDir filesep 'paramsIn.mat'],'p');
    
    display(['Finished Neurite Orientation Consistency Test : ' movieData.outputDirectory_ ' for Ch' num2str(iCh) ]); 
end % iCh
end % The END

