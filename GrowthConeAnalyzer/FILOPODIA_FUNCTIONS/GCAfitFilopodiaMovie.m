function [ filoInfo] = GCAfitFilopodiaMovie(movieData,paramsIn)
%fitLinescansMovie: performs automated fitting of filopodia detections

% INPUT:
% filoInfo: output of filopodiaReconstruct function
%           the N filo long structure- an output of the filopodia reconstruct
%           provides both steerable filter based reconstructions and the forward
%           projections used for fitting
%
% paramsIn -Structure with inputs for optional parameters. The
%           parameters should be stored as fields in the structure, with the field
%           names and possible values as described below
%
% Input/output Generic Wrapper Params:
%       ('OutputDirectory' -> character string) Optional. A character
%       string specifying the directory to save the analInfo with the new
%       filopodia fits
%       If not input, the vectors will be saved to the same directory as
%       the movieData, in a sub-directory called "filopodia_fits"
%
%       ('ChannelIndex' -> Positive integer scalar or vector) Optional. The
%       integer index of the channel(s) to fit the filopodia. If not input,
%       all channels with a filopodia reconstruction will be used.
%
% Internal Function Params:
% ('InternalFiloOn' -> scalar ) Optional If true , fit internal filo as marked in filoInfo data struct,
%      if false, fits external filo as marked in data struct
%      NOTE TO SELF: (need to have an option where it will run through both!!!)
%
% ('NumPixForFitBack' -> scalar) Optional Default = 10 Pixels
%         This parameter dictates the number of pixels back along the filopodia
%         that will be used in the fit relative to the end of the high confidence filopodia tip
%         estimated via the steerable filter response thresholding.
%         Note the signal is often quite noisy along the filopodia - ie there can be
%         multiple possible viable sigmoidal fits. This is especially the
%         case when fitting a lifeAct reporter signal- so one just wants to fit the
%         signal in a local area around the putative tip of the filopodia
%         If the number of pixels in the filopodia is less than this value
%         the entire filopodia length from the thresholded steerable filter response
%         will be used for the respective filopodia tip localization fit.
%
% NOTE TO SELF: this was actually dictated first in the original walkFiloForandBack function
%         ('NumberPixForFitFor' -> scalar) Optional Default = 10 Pixels
%         This parameter dictates the number of pixels forward (relative to
%         the steerable filter threshold set in the filopodia
%         reconstruction) used in the fitting.
%
%
%
% ('ValuesForFit' -> character 'Intensity','Response','Both') Optional
% Default = 'Intensity'
%         The values to use the fitting:
%         Intensity: The values corresponding to the intensity will be fit
%         to a sigmoid
%         Response: The values corresponding to the NMS response from the
%         steerable filter will be fit to a sigmoid.
%         Both: Both of the above operations will be performed- this is
%         mainly for comparison in the early stages of development
%
%('Averaging' -> character ('Gaussian Weighted' , 'Perpendicular'
% For Fitting: Response or
%

%
%('SavePlots' -> logical) Optional Default = 1
%
%
%
%
%
%


% Output:
%
%   movieData - the updated MovieData object with the
%   parameters, paths etc. stored in it, in the field movieData.processes_.
%
%   The analInfo is written to the directory specified by the
%   parameter OutputDirectory, they are stored as a single .mat file.
%   filoInfo will be updated with fields...
%
%% ----- Input ------ %%
if nargin < 1 || ~isa(movieData,'MovieData')
    error('The first input must be a valid MovieData object!')
end

if nargin < 2
    outputDirectory = [movieData.outputDirectory_ filesep 'filopodia_fits'];
    paramsIn.OutputDirectory = outputDirectory ;
    paramsIn.ChannelIndex = 1  ;
    paramsIn.InternalFiloOn = 3; % types 1, 2, or 3.
    paramsIn.NumPixForFitBack = 10; % should maybe eventually make this distance?
    paramsIn.ValuesForFit = 'Intensity'; % default is the intensity;
    paramsIn.SavePlots = 1;
    paramsIn.restart.startFrame = 'auto';   % 'auto'
    paramsIn.restart.endFrame = 'auto'; % 'auto' or number, default auto.
end
p = paramsIn ;
%Get the indices of any previous mask refinement processes from this function
% iProc = movieData.getProcessIndex('FitFiloProcess',1,0);
%
% %If the process doesn't exist, create it
% if isempty(iProc)
%     iProc = numel(movieData.processes_)+1;
%     movieData.addProcess(ProtrusionProcess(movieData,movieData.outputDirectory_));
% end
%
% %Parse input, store in parameter structure
% p = parseProcessParams(movieData.processes_{iProc},paramsIn);
%% Start Wrapper
for iCh = 1:numel(paramsIn.ChannelIndex)
    
    % Make Output Directory
    outPutDirC = [p.OutputDirectory filesep 'Filopodia_Fits_Channel_' num2str(p.ChannelIndex(iCh)) ];
    %% Get Start and End Frame Based on Restart Choices
    withFiloFit  = [outPutDirC filesep 'analInfoTestSave.mat'];
    % If file exists
    if  exist(withFiloFit,'file')==2;
        load(withFiloFit) % load the file
        display('Loading Previously Run Filopodia Fits');
        if strcmpi(p.restart.startFrame,'auto')
            startFrame = find(arrayfun(@(x) ~isfield(analInfo(x).filoInfo, 'Ext_exitFlag')...
                ,1:length(analInfo)),1,'first');
            startFrame  = startFrame-1;
            display(['Auto Start: Starting Neurite Body Reconstruction at Frame ' num2str(startFrame)]);
        else
            startFrame = p.restart.startFrame; % use user input
            display(['Manual Start: Starting Neurite Body Reconstruction at Frame ' num2str(startFrame)]);
        end
    else % if doesn't exist : load from reconstruct file
        withoutFiloFit =  [movieData.outputDirectory_ filesep 'filopodia_reconstruct'  filesep ...
            'Filopodia_Reconstruct_Channel_' num2str(p.ChannelIndex(iCh)) filesep 'analInfoTestSave.mat'];
        % load reconstruction data
        load(withoutFiloFit) ; %
        %  test to make sure
        firstEmpty= find(arrayfun(@(x) isempty(x.reconstructInfo),analInfo),1,'first');
        if firstEmpty ~= movieData.nFrames_;
            
            allEmpty = find(arrayfun(@(x)  isempty(x.reconstructInfo),analInfo));
            display('Note: Filopodia reconstructions have not been performed for all frames of movie')
            
            arrayfun(@(x) display(['No Filopodia Reconstruction Found for Frame ' num2str(allEmpty(x))]),1:length(allEmpty));
        end % firstEmpty
        startFrame = 1; 
    end % if exist withFiloFit
    
    if strcmpi(p.restart.endFrame,'auto');
        endFrame = movieData.nFrames_;
        display(['Auto End: Fitting Filopodia From Frame ' num2str(startFrame) ' to ' num2str(endFrame)]);
    else
        endFrame = p.restart.endFrame;
        display(['Manual End: Fitting Filopodia From Frame ' num2str(startFrame) ' to ' num2str(endFrame)]);
    end
    
    %% Start Loop Over Movie
    % GET FRAME INFORMATION - this function wraps per frame
    for iFrame = startFrame:endFrame
        % get the filoInfo for the current frame
        filoInfo = analInfo(iFrame).filoInfo;
        imgPath = [movieData.getChannelPaths{p.ChannelIndex(iCh)} filesep movieData.getImageFileNames{p.ChannelIndex(iCh)}{iFrame}];
        img = double(imread(imgPath));
        % make a specific output directory for the plotting for each frame
        pSpecific = p;
        pSpecific.sigma = movieData.channels_.psfSigma_;
        if isempty(pSpecific.sigma)
            display(['Using sigma 0.43']);
            pSpecific.sigma = 0.43;
        end
        if pSpecific.SavePlots == 1
            pSpecific.OutputDirectory = [outPutDirC filesep 'Linescans' filesep 'Frame ' num2str(iFrame,'%03d')];
            mkClrDir(pSpecific.OutputDirectory)
        end
        
        filoInfo = GCAfitFilopodia(filoInfo,img,pSpecific) ;
        % rewrite the filoInfo with the extra filo Info fields.
        analInfo(iFrame).filoInfo = filoInfo;
        display(['Finished Fitting Filopodia for  Channel ' num2str(p.ChannelIndex(iCh)) 'Frame ' num2str(iFrame)]);
        analInfo(iFrame).reconstructInfo.createTimeFiloFit = clock;
        hashTag = gcaArchiveGetGitHashTag;
        analInfo(iFrame).reconstructInfo.hashTagFiloFit = hashTag;
        
        save([outPutDirC filesep 'analInfoTestSave.mat'],'analInfo','-v7.3')
        
    end % for iFrame
end % for iCh
