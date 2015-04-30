function [analInfo] = GCAReconstructFilopodiaMovie(movieData,paramsIn)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUT
%
%   movieData - The MovieData object describing the movie, as created using
%   movieSelectorGUI.m
%
%   paramsIn - Structure with inputs for optional parameters. The
%   parameters should be stored as fields in the structure, with the field
%   names and possible values as described below:
%
%
%   Parameter Structure Field Names:
%
% Generic Fields: (Input/Output Fields Needed for Wrapper)
%       ('OutputDirectory' -> character string) Optional. A character
%       string specifying the directory to save the analInfo structure to.
%       If not input, the analInfo will be saved in the same directory
%       as the movieData, in a sub-directory called
%       "neurite_orientation_estimations"
%       % PERSONAL NOTE : You might need
%       to truncate these names for the windows paths. %
%
%       ('ChannelIndex'-> Positive integer scalar or vector) The integer
%       indices of the channel(s) on which to perform the neurite orientation
%       estimation.
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
%  GCAReconstructFilopodia specific functions
%
% OUTPUT: (see main function GCAReconstructFilopodia- for details):
%       analInfo Rx1 structure array where R = the number of frames
%       with added field .filoInfo a Rx1 structure where R is the number of
%       filopodia detected in that frame 
%       .filoInfo has many subfields that store information regarding the
%       filopodias coordinates,orientation,grouping(for branching) etc. 

p = paramsIn;



%% Init:
nFrames = movieData.nFrames_;
nChan = numel(p.ChannelIndex);
imSize = movieData.imSize_;
ny = imSize(1);
nx = imSize(2);

%%
%% Loop for each channel
for iCh = 1:nChan
    
    display(['Finding Neurite Orientation Channel ' num2str(iCh)]);
       
    %% Get Start and End Frames Based on Restart Choice
    
    
    
    outDir = [movieData.outputDirectory_ filesep ...
        'filopodia_reconstruct' filefsep 'Filopodia_Reconstruct_Channel_' num2str(iCh)]; 
    
    if ~isdir(outDir);
        mkdir(outDir);
    end 
    
    reconstructFile = [outDir filesep 'analInfoTestSave.mat']; 
      
    % If file exists
    if  exist(reconstructFile,'file')==2;
        load(reconstructFile) % load the file
        display('Loading Previously Run Filopodia Reconstructions ');
        if strcmpi(p.restart.startFrame,'auto')
            startFrame = numel(analInfo)-1;
            if startFrame == 0
                startFrame = 1; % reset to 1;
            end
            display(['Auto Start: Starting Filopodia Reconstruction at Frame ' num2str(startFrame)]);
        else
            startFrame = p.restart.startFrame; % use user input
            display(['Manual Start: Starting Filopodia Reconstruction at Frame ' num2str(startFrame)]);
        end
    else % if doesn't exist
        startFrame = 1;
        display('No Filopodia Reconstruction Folder Found: Creating and Starting at Frame 1');
        % load analInfo from neuriteEstimation folder 
        load( [movieData.outputDirectory_ filesep 'neurite_body_masks' filesep 'Neurite_Body_Masks_Channel_1' ... 
            filesep 'analInfoTestSave.mat']);    
    end % exist(orientFile,'file') == 2
    
    if strcmpi(p.restart.endFrame,'auto');
        endFrame = nFrames;
        display(['Auto End: Performing Filopodia Reconstructions From Frame ' num2str(startFrame) ' to ' num2str(endFrame)]);
    else
        endFrame = p.restart.endFrame;
        display(['Manual End: Performing Filopodia Reconstructions From Frame ' num2str(startFrame) ' to ' num2str(endFrame)]);
    end
    %%  
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
     
%% Test if the protrusion process was run 
     idxProt = find(cellfun(@(x) sum(strcmpi(x.name_,'Protrusion')),MD.processes_));
    if ~isempty(idxProt)
        % load the three protrusion process associated cell structures, 
        % 'normals','protrusion','smoothedEdge'
        load(MD.processes_{idxProt(end)}.outFilePaths_);
        % If for some reason there is more than one protrusion process 
        % associated with the data, tell the user that your are using the
        % most recently run
        if length(protS) >1 
            display(['There was more than one veil protrusion process associated with ' ... 
                'this movie: local filopodia to veil calculations will be '... 
                'will be performed using the most recently run protrusion process']); 
        end % if length(protS)           
    else % no veil protrusion vectors were found
        display(['The veil protrusion process was not found for this movie :'... 
            'The filopodia reconstruction will continue, but no local filopodia orientation ' ... 
            'calculations relative to the veil will be performed']); 
        normals = []; 
        smoothedEdge = []; 
    end % ~isempty(idxProt)
     
    
    
    
    %% Start Movie Loop %%%%
    for iFrame = startFrame:endFrame
        % Load image
        img = double(imread( [listOfImages{iFrame,2} filesep listOfImages{iFrame,1}] ));
        
        if ~isempty(normals)
        protrusionC.normal = normals{iFrame}; 
        protrusionC.smoothedEdge = smoothedEdge{iFrame}; 
        else 
            protrusionC = []; 
        end 
        
        [analInfo,troubleshootFigs] = GCAReconstructFilopodia(img,analInfoC,protrusionC,paramsIn);
        
        
        
        
        
        % quick fix for the plots is to just make the frame number an input for not 20140812
        hashTag =  gcaArchiveGetGitHashTag;
        analInfo.hashTag = hashTag; % make sure to add the hash tag first so the structure is similar (or initiate in begin)
        analInfo(iFrame) = analInfo;
        save( [outDir filesep 'analInfo.mat'],'analInfo');
    end % iFrame
    
end   % for iCh
end


