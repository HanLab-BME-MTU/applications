function [ output_args ] = GCAlabelFiloAccordingToVeilBehaviorMovie(movieData,paramsIn)
%Movie Wrapper for GCAlabelFiloAccordingToVeilBehavior (start with just the
% persistence option
% 
%% CHECK Parameters 
if nargin < 1 || ~isa(movieData,'MovieData')
    error('The first input must be a valid MovieData object!')
end

if nargin < 2
    % Generic
    paramsIn.OutputDirectory = [movieData.outputDirectory_ filesep 'EXPLORE_VEIL_FILO_INTERACT'];
    paramsIn.ChannelIndex = 1;
    paramsIn.ProcessIndex = 0; % use raw images
    % for now paramsIn.GroupType will allow one to turn on different
    % labels that will add fields to the filoData struct
    paramsIn.GroupType.Persistence.flag = 1; 
    
    
    paramsIn.GroupType.Persistence.cutoff = 120; % in sec 
    
    paramsIn.labelMovie = 1; 
end 

p = paramsIn;

%% Init:
% make the new directory to begin to explore the spatial patterns 
mkClrDir(p.OutputDirectory); 
if p.GroupType.Persistence.flag ==1 
    saveDir  = [p.OutputDirectory filesep 'PersistentVeil']; 
    mkClrDir(saveDir);  %maybe change later. 
    
end 
% if the movie feature is selected
 if p.labelMovie ==1 
     % if the persistence feature is seleteced 
     if p.GroupType.Persistence.flag == 1
     movieDir = [saveDir filesep 'Movie']; 
     mkClrDir(movieDir)
     end 
 end 
 
     
    
nFrames = movieData.nFrames_;
nChan = numel(p.ChannelIndex);
imSize = movieData.imSize_;
ny = imSize(1); 
nx = imSize(2); 
% eventually make protrusion and retraction different 
if p.GroupType.Persistence.flag == 1 % eventually make protrusion and retraction
    
    
    % load edge velocity data
    load([ movieData.outputDirectory_ filesep 'edgeVelocityAnalysis' filesep 'EdgeMotion.mat']);
    
    output = getWindowsWithHighPersistence(analysisResults,p.GroupType.Persistence.cutoff); % time series so whole movie
   
    
    
end


for iCh = 1:nChan
  if p.labelMovie == 1 
      imgNames = movieData.getImageFileNames(iCh); 
  end    
    
  % INSERT CHECKS TO MAKE SURE THE FILO RECONSTRUCT WAS RUN  
  % load analInfo this contains the filoInfo for each frame    
  load([movieData.outputDirectory_ filesep 'filopodia_reconstruct' filesep... 
      'Filopodia_Reconstruct_Channel_' num2str(iCh)' filesep 'analInfoWithWindInfo.mat']); 
  display(['Labeling Filopodia Surrounded By High Persistence Veil Motion Channel ' num2str(iCh)]); 
 % INSERT CHECKS TO MAKE SURE THE WINDOWING WAS RUN- Likely a nice way to
 % do this with  movieData to get all files associated with the windowing
 % step- for now
 for iFrame = 1:nFrames-1
     
    if p.labelMovie ==1 
        img = double(imread([ movieData.channels_(iCh).channelPath_ filesep imgNames{1}{iFrame}]));
        % load the windows 
        load([movieData.outputDirectory_ filesep 'windows' filesep 'windows_frame__frame_' num2str(iFrame,'%02d') '.mat']); 
        setFigure(nx,ny,'on'); 
        imshow(-img,[]); 
        hold on 
    end 
     
    
     filoInfo = analInfo(iFrame).filoInfo; 
     
     filoInfo = markFiloWithPersistentVeil(output,filoInfo,iFrame);
      %% make a small matrix of the windows marked by persistent retraction
      % want to make an idx matrix so can mark which ones have high persistence 
     % for each window test if pers retract if yes 1 if no 0 
     %nWinds = length(output.retractionAnalysis.windows); 
   
    % frameNumPers= output.(fieldC).windows(filoInfo(idxC).windowIdx).blockOut; 
    
    % extract from cells (cell number tells block idx - we don't care at
    % this point) 
   % frameNumPers = vertcat(frameNumPers{:});
   % currentFrame= sum(frameNumPers==iFrame) ;
     
   %%
     if p.labelMovie ==1 
         % plot the windows and the filoAssociated with persistent
         % protrusion/retraction   
         persType{1} = 'protrusion' ;
         persType{2} = 'retraction'; 
         c = ['g','c']; 
         
         for i = 1:2 
         filoPlotC = filoInfo(logical(vertcat(filoInfo(:).([persType{i} 'persVeil']))));
             
         plotfilosIntAndExt(filoPlotC,[ny,nx],1,1,c(i),0);     
         
         end
          plotWindows(windows,{'g','FaceAlpha',0},'bandMax',1)
          nWinds = length(output.retractionAnalysis.windows);
          idxPerRetract =  arrayfun(@(x) sum(vertcat(output.retractionAnalysis.windows(x).blockOut{:})==iFrame),1:nWinds); 
          plotWindows(windows(logical(idxPerRetract)),{'b','FaceAlpha',0.5},'bandMax',1); 
          
          idxPerProt  = arrayfun(@(x) sum(vertcat(output.protrusionAnalysis.windows(x).blockOut{:})==iFrame),1:nWinds); 
          
          plotWindows(windows(logical(idxPerProt)),{'r','FaceAlpha',0.5},'bandMax',1); 
          
        roiYX=  bwboundaries(analInfo(iFrame).masks.neuriteEdge);
        cellfun(@(x) plot(x(:,2),x(:,1),'Color','y'),roiYX);  
         saveas(gcf,[movieDir filesep num2str(iFrame, '%03d') '.png']);
     saveas(gcf,[movieDir filesep num2str(iFrame, '%03d') '.fig']); 
     close gcf
     end 
         
  %reWrite the analInfo 
analInfo(iFrame) = filoInfo;    
          
     
 end % iFrame


save([saveDir filesep 'analInfo.mat']); 




end
end 
