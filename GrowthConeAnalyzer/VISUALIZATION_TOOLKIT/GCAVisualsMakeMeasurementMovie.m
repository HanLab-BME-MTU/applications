function [ output_args ] = GCAVisualsMakeMeasurmentMovie(MD,varargin)
% GCAVisualsMakeMeasurementMovie 
% 
% 
% This function will make the measurement movies for all measurements 
% Default is to ask the user which type of movie they would like to make 

% INPUT: MD movieData object
% 
%% Check input
ip = inputParser;

ip.CaseSensitive = false;

ip.addParameter('interactive',true,@(x) islogical(x));

% Note need to make more generic for channel wrap! 
defaultMeasDir = [MD.outputDirectory_ filesep 'MEASUREMENT_EXTRACTION' ]; 
ip.addParameter('MeasurementDirectory',defaultMeasDir,@(x) ischar(x)); 

defaultInDir = [MD.outputDirectory_ filesep 'SegmentationPackage' ... 
    filesep 'StepsToReconstruct' filesep 'VII_filopodiaBranch_fits' filesep 'Channel_1']; 

defaultVeilDir = [MD.outputDirectory_ filesep 'SegmentationPackage' ... 
    filesep 'StepsToReconstruct' filesep 'IV_veilStem_length' filesep 'Channel_1']; 


ip.addParameter('InputDirectory', defaultInDir,@(x) ischar(x)); 

ip.addParameter('ScaleBar',true,@(x) islogical(x)); 
ip.addParameter('TreatmentFrame',[]); 
ip.addParameter('TreatmentTitle','CK666 Added'); 
ip.addParameter('veilDir',defaultVeilDir); 
ip.addParameter('plotText',true); 

ip.parse(varargin{:});
p = ip.Results;

%% go into each folder and look for a measurement 
% for now just search files - redesign so that the parameters in the
% future are more cleverly named

% might also include ylabel name and ylim for each parameter and
% read in each Descriptor directory to keep constant.

% search all descriptor parameters.
localParamFiles = searchFiles('meas_',[],[ip.Results.MeasurementDirectory filesep 'Descriptor'],1);

paramNamesC = cellfun(@(x) strrep(x,'meas_',''),localParamFiles(:,1),'uniformoutput',0);
paramNamesC = cellfun(@(x) strrep(x,'.mat',''),paramNamesC,'uniformoutput',0);

% if interactive
% ask the user for which measurement they would like to make a sanity movie

% most measurments are going to be simple it will just be plotting
% the filopodia of the filter set with the value saved
if ip.Results.interactive == true
    paramSelect  = listSelectGUI(paramNamesC,[],'move');
    selected  = paramNamesC(paramSelect);
else % make all movies
    selected = paramNamesC; 
end

% load the filoInfo
load([ip.Results.InputDirectory filesep 'filoBranch.mat']);  
imgSize = MD.imSize_; 

  % plot the veil in black 
       load([ip.Results.veilDir filesep 'veilStem.mat']); 

%% Different File types might require different functions
for iSelect = 1:numel(selected)
    
    % make the directory to save the measurement movie 
    outDir = [localParamFiles{paramSelect(iSelect),2} '_' selected{iSelect} '_Measurement_Movie']; 
    if ~isdir(outDir) 
        mkdir(outDir) ; 
    end 
    
    
    load([localParamFiles{paramSelect(iSelect),2} filesep localParamFiles{paramSelect(iSelect),1}]);
    x = upDirectory(localParamFiles{paramSelect(iSelect),2},1);
    load([x filesep 'filoFilterSet.mat']);
    
    
    
    for iFrame = 1:MD.nFrames_-1
        filoInfo = filoBranch(iFrame).filoInfo;
        plotValues = measC{iFrame};
        filterSetC= filoFilterSet{iFrame};
        img = double(imread([MD.getChannelPaths{1} filesep MD.getImageFileNames{1}{iFrame}]));
        
        
        setFigure(imgSize(2),imgSize(1),'on');
        
        imshow(-img,[]) ;
        hold on
        
        %testForBranchMode = (~isempty(regexpi(selected{iSelect},'branchDistanceTo'))  || ~isempty(selected{iSelect})); 
        % testForBranchMode = (~isempty(regexpi(selected{iSelect},'branchDistanceTo'))  || ~isempty(selected{iSelect})); 
        testForBranchMode = false; 
        if testForBranchMode == 1
            branchMode = true;
        else
            branchMode = false;
        end
       
      % Note: maybe try to fix this input before release - historical and a
      % bit cumbersome. 
      if ~isempty(regexpi(selected{iSelect},'Embedded')); 
          filoPlotType  = 2; 
      elseif ~isempty(regexpi(selected{iSelect},'FullActinBundle'));
          filoPlotType = 3; 
      else 
          filoPlotType = 1; 
      end 
          
    
        
       GCAVisualsFilopodiaMeasurementOverlays(filoInfo,imgSize,... 
            'filoFilterSet',filterSetC,'plotValues',plotValues,...
            'branchMode',branchMode,'colorByValue',true,'plotText',ip.Results.plotText,'justExt',filoPlotType);
       hold on 
       if ip.Results.ScaleBar == true 
       pixSizeMic = MD.pixelSize_./1000; 
       plotScaleBar(10/pixSizeMic,5,'Color',[0,0,0],'Location','SouthEast');
       end 
        
       if ~isempty(ip.Results.TreatmentFrame) 
           if iFrame >=ip.Results.TreatmentFrame; 
               text(10,10,ip.Results.TreatmentTitle,'color','k'); 
           end
       end 
       
     
       veilStemMask = veilStem(iFrame).finalMask; 
       roiYX = bwboundaries(veilStemMask); 
       cellfun(@(x) plot(x(:,2),x(:,1),'color','k'),roiYX);  
       
       saveas(gcf,[outDir filesep num2str(iFrame,'%03d') '.png']); 
       
       close gcf 
        
    end


end

