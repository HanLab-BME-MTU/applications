function [ filoInfo] = GCAfitFilopodiaMovie( movieData,paramsIn)
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
% Default = 'Both'
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
    paramsIn.InternalFiloOn = 0; 
    paramsIn.NumPixForFitBack = 10; % should maybe eventually make this distance?
    paramsIn.ValuesForFit = 'Intensity'; % default is the intensity;
    paramsIn.SavePlots = 0; 
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
   
        outPutDirC = [p.OutputDirectory filesep 'Filopodia_Fits_Channel_' num2str(p.ChannelIndex(iCh))]; 
    
        
        mkClrDir(outPutDirC) 
        
        filoInfoDir  = [movieData.outputDirectory_ filesep  'filopodia_reconstruct' filesep 'Filopodia_Reconstruct_Channel_' num2str(iCh)];
     % load reconstruction data including filoInfo per frame 
    load([filoInfoDir filesep 'analInfoTestSave.mat']) ; %
    

    
    
  % GET FRAME INFORMATION - this function wraps per frame 
    for iFrame = 1:length(analInfo)
        % get the filoInfo for the current frame
        filoInfo = analInfo(iFrame).filoInfo;
        imgPath = [movieData.getChannelPaths{p.ChannelIndex(iCh)} filesep movieData.getImageFileNames{p.ChannelIndex(iCh)}{iFrame}];
        img = double(imread(imgPath)); 
        % make a specific output directory for the plotting for each frame 
        pSpecific = p; 
        pSpecific.sigma = movieData.channels_.psfSigma_; 
        if isempty(pSpecific.sigma) 
            display(['Using sigma 0.5']); 
            pSpecific.sigma = 0.5;
        end 
        if pSpecific.SavePlots == 1 
            
        pSpecific.OutputDirectory = [outPutDirC filesep 'Linescans' filesep 'Frame ' num2str(iFrame,'%03d')];
        mkClrDir(pSpecific.OutputDirectory)
        end 
        
        filoInfo = GCAfitFilopodia(filoInfo,img,pSpecific) ; 
        % rewrite the filoInfo with the extra filo Info fields. 
        analInfo(iFrame).filoInfo = filoInfo;
        display(['Finished Fitting Filopodia for  Channel ' num2str(p.ChannelIndex(iCh)) 'Frame ' num2str(iFrame)]); 
       
        save([outPutDirC filesep 'analInfoTestSave.mat'],'analInfo','-v7.3')
        
        
     
        
        

        
       
            
            
         
           
            
                
                
                
                
                
              
                
                
               
                    
                    
                    
                    
                
                
                
                
                
              
           
 
    end % for iFrame
end % for iCh










