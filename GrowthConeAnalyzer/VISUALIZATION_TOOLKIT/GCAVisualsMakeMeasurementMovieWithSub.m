function [ output_args ] = GCAVisualsMakeMeasurementMovie(MD,varargin)
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
ip.addParameter('measurements',[],@(x) iscell(x));

defaultFrames= 1:MD.nFrames_ -1;
ip.addParameter('frames', defaultFrames,@(x) isnumeric(x));


ip.addParameter('MeasurementDirectory',[],@(x) ischar(x) || isempty(x));


ip.addParameter('InputDirectory', [],@(x) ischar(x) || isempty(x));

ip.addParameter('ScaleBar',false,@(x) islogical(x) );

ip.addParameter('TreatmentFrame',[]);
ip.addParameter('TreatmentTitle','CK666 Added');

ip.addParameter('VeilDirectory',[]);
ip.addParameter('plotText',false);


ip.addParameter('cMapLimits',[]); 
ip.addParameter('colorByValue',true);
ip.addParameter('extraColor',[1 1 1]); % adds an extra color to the color bar so that NaN values are white

ip.addParameter('visible','off'); % if on show figure while plotting
ip.addParameter('colorbarOverlay',true);

ip.addParameter('otherFiles',true);
ip.addParameter('SubRegionFlag',false,@(x) islogical(x))


defaults{1,1} = 'filoDensityAlongVeil'; defaults{1,2} = [0,10];
defaults{2,1} = 'filoOrient'; defaults{2,2} = [0,180];
defaults{3,1} = 'filoIntensityEmbedded';defaults{3,2} = [0.5,2];
defaults{4,1} = 'filoIntensityToVeil'; defaults{4,2} = [0.5,2];
defaults{5,1} = 'filoLengthEmbedded';defaults{5,2} = [0,10];
defaults{6,1} = 'filoLengthFullActinBundle';defaults{6,2} = [0,10];
defaults{7,1} = 'filoLengthToVeil'; defaults{7,2} = [0,10];
defaults{8,1} = 'filoCurvature'; defaults{8,2} = [0,.5];
defaults{9,1} = 'branchLength_2ndOrder'; defaults{9,2} =[0,10];
defaults{10,1} = 'branchOrientation_2ndOrder' ; defaults{10,2} = [0,180];
defaults{11,1} = 'branchIntensity_2ndOrder' ; defaults{11,2} = [0 2];
defaults{12,1} = 'validation' ; defaults{12,2} = [0,10];

ip.addParameter('minMaxDefaults',defaults); % defaults for me are set below

ip.parse(varargin{:});
p = ip.Results;
%% 
if isempty(ip.Results.MeasurementDirectory)
    if ~ip.Results.SubRegionFlag
        % Note need to make more generic for channel wrap!
        measDir = [MD.outputDirectory_ filesep 'MEASUREMENT_EXTRACTION' ];
    else
        measDir = [MD.outputDirectory_ filesep 'GCAMeasurementExtraction' filesep 'GC'];
    end
else
    measDir = ip.Results.MeasurementDirectory;
end
    
if isempty(ip.Results.InputDirectory)
    if ~ip.Results.SubRegionFlag
        inDir = [MD.outputDirectory_ filesep 'SegmentationPackage' ...
            filesep 'StepsToReconstruct' filesep 'VII_filopodiaBranch_fits' filesep 'Channel_1'];
    else
        inDir = [MD.outputDirectory_ filesep 'SegmentationPackage' ...
            filesep 'GCASubRegions' filesep 'GC' filesep 'filopodia'];
    end
else
    inDir = ip.Results.InputDirectory;
end


if isempty(ip.Results.VeilDirectory)
    if ~ip.Results.SubRegionFlag
        veilDir = [MD.outputDirectory_ filesep 'SegmentationPackage' ...
            filesep 'StepsToReconstruct' filesep 'IV_veilStem_length' filesep 'Channel_1'];
    else
        veilDir = [MD.outputDirectory_ filesep 'SegmentationPackage' filesep 'GCASubRegions' filesep ...
            'GC' filesep 'masks'];
    end
else 
    veilDir = ip.Results.VeilDirectory;
end

%% Set up my personal default range for each parameter set as this will keep
% everything constant for comparisons and ensure I can simply loop through
% all the visualization overlay for any given movie

%% go into each folder and look for a measurement
% for now just search files - redesign so that the parameters in the
% future are more cleverly named

% might also include ylabel name and ylim for each parameter and
% read in each Descriptor directory to keep constant.

% search all descriptor parameters.
localParamFiles = searchFiles('meas_',[],[measDir filesep 'Descriptor'],1);

paramNamesC = cellfun(@(x) strrep(x,'meas_',''),localParamFiles(:,1),'uniformoutput',0);
paramNamesC = cellfun(@(x) strrep(x,'.mat',''),paramNamesC,'uniformoutput',0);

%% ACTION REQUIRED 
       % NOTE : as of 20160125 need to add the ID plotting option 
%         paramNamesC{end+1,1} = 'IDs';
%         localParamFiles{end+1,2} = [measDir filesep 'Descriptor' filesep 'IDs'];
%% 

% if interactive
% ask the user for which measurement they would like to make a sanity movie

% most measurments are going to be simple it will just be plotting
% the filopodia of the filter set with the value saved
if ip.Results.interactive == true
  
    paramSelect  = listSelectGUI(paramNamesC,[],'move');
    selected  = paramNamesC(paramSelect);
      
else % make all movies
    if ~isempty(ip.Results.measurements);
 
             
        selected = ip.Results.measurements;
        % find those measurements from the parameter names
        test = arrayfun(@(x) strcmpi(selected{x},paramNamesC),1:numel(selected),'uniformoutput',0);
        paramSelect = cellfun(@(x) find(x),test);
       
        %paramSelect  = find(sum(horzcat(test{:}),2));
    else % loop through all found
        selected = paramNamesC ;
        paramSelect = 1:numel(selected);
    end
end






% load the filoInfo
load([inDir filesep 'filoBranch.mat']);
imgSize = MD.imSize_;

if ~ip.Results.SubRegionFlag
    % plot the veil in black
    load([veilDir filesep 'veilStem.mat']);
else % note need to fix this output to remain consistent either one or the other formate
    x = getFileListFromFolder(veilDir);
    
    for iFrame = 1:numel(x)
        mask = logical(imread(x{iFrame})); 
        veilStem(iFrame).finalMask = mask ;
        clear mask
    end
    
    if sum(cellfun(@(x) strcmpi(x,'filoOrient'),selected));
        % load from default veil stem directory %% NOTE
        % maybe change design of the subRegional mask
        % output to accomodate this later.
        x = load([MD.outputDirectory_ filesep 'SegmentationPackage' filesep...
            'StepsToReconstruct' filesep 'IV_veilStem_length' filesep ...
            'Channel_1' filesep 'veilStem.mat']);
        for iFrame = 1:length(x.veilStem);
            backbone = x.veilStem(iFrame).neuriteLongPathIndices;
            veilStem(iFrame).neuriteLongPathIndices = backbone;
        end
    end
end
    
%%
nFrames = length(ip.Results.frames);
frames = ip.Results.frames;
% check to make sure the filter and the number of frames is equivalent

%% Different File types might require different functions
for iSelect = 1:numel(selected)
  if strcmpi(paramSelect,'IDs')
      if ~isdir(localParamFiles{paramSelect(iSelect),2}); 
          mkdir(localParamFiles{paramSelect(iSelect),2}); % make the ID directory 
      end
  end 
  
 
    % make the directory to save the measurement movie
    outDir = [localParamFiles{paramSelect(iSelect),2} filesep selected{iSelect} '_Measurement_Movie'];
    if ~isdir(outDir)
        mkdir(outDir) ;
    end
    
   
    load([localParamFiles{paramSelect(iSelect),2} filesep localParamFiles{paramSelect(iSelect),1}]);
    
    x = upDirectory(localParamFiles{paramSelect(iSelect),2},1);
    
    load([x filesep 'filoFilterSet.mat']);
    
    
    for iFrame = 1:nFrames
        frameC = frames(iFrame);
       
            filoInfo = filoBranch(frameC).filoInfo;
            
            plotValues = measC{frameC};
            if ~isempty(plotValues);
                
                if strcmpi(ip.Results.cMapLimits,'defaults')
                    [cMapLimits] =  defaults{strcmpi(selected{iSelect},defaults(:,1)),2} ;
                elseif isempty(ip.Results.cMapLimits)
                    cMapLimits(1) = min(plotValues);
                    cMapLimits(2) = max(plotValues);
                end
                
                
                
                filterSetC= filoFilterSet{frameC};
                %
                if ~isempty(regexpi(selected{iSelect},'Embedded'));
                    filterSetC = (filterSetC(:,1)==1 & filterSetC(:,2) ==1);
                end
                img = double(imread([MD.getChannelPaths{1} filesep MD.getImageFileNames{1}{frameC}]));
                
                
                setFigure(imgSize(2),imgSize(1),ip.Results.visible);
                
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
                
                
                
                % Know this is a bit cumbersome but for now let's just have a flag to
                % redirect if need to plot curvature
                if ~strcmpi(selected{iSelect},'filoCurvature')
                    
                    
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
                        'branchMode',branchMode,'colorByValue',ip.Results.colorByValue,'plotText',ip.Results.plotText,'justExt',...
                        filoPlotType,'extraColor',ip.Results.extraColor,'cMapLimits',cMapLimits);
                    
                    if strcmpi(selected{iSelect},'filoOrient') ;
                            overlay = zeros(imgSize);
                            backbone = veilStem(frameC).neuriteLongPathIndices;
                            overlay(backbone) = 1;
                            spy(overlay,'k');
                    end
                    
                    
                else
                    
                    GCAVisualsColorCodeByCurvature(filoInfo,'filoFilterSet',filterSetC,'cMapLimits',cMapLimits);
                    
                end
                
                hold on
                if ip.Results.ScaleBar == true
                    pixSizeMic = MD.pixelSize_./1000;
                    plotScaleBar(10/pixSizeMic,2,'Color',[0,0,0],'Location','SouthEast');
                end
                
                if ~isempty(ip.Results.TreatmentFrame)
                    if frameC >=ip.Results.TreatmentFrame;
                        text(10,10,ip.Results.TreatmentTitle,'color','k');
                    end
                end
                
               
                veilStemMask = veilStem(frameC).finalMask;
                roiYX = bwboundaries(veilStemMask);
                cellfun(@(x) plot(x(:,2),x(:,1),'color','k'),roiYX);
                
                saveas(gcf,[outDir filesep num2str(frameC,'%03d') '.png']);
                
                if ip.Results.otherFiles
                    saveas(gcf,[outDir filesep num2str(frameC,'%03d') '.eps'],'psc2');
                    saveas(gcf,[outDir filesep num2str(frameC,'%03d') '.fig']);
                end
                
                close gcf
                
                
                if ip.Results.colorbarOverlay
                    figure('visible',ip.Results.visible)
                    test = zeros(imgSize);
                    test(1,1) = cMapLimits(1);
                    test(1,2) = cMapLimits(2);
                    imagesc(test);
                    colormap(jet(128));
                    colorbar
                    
                    saveas(gcf,[outDir filesep 'ColorBarOverlay_' selected{iSelect} '.eps'],'psc2');
                    saveas(gcf,[outDir filesep 'ColorBarOverlay_' selected{iSelect} '.fig']);
                    close gcf
                    
                end
            end
    end 
    
end

