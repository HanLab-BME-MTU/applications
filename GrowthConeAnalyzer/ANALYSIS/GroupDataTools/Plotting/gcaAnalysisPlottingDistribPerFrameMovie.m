function [h] = gcaAnalysisPlottingDistribPerFrameMovie(MD,varargin)
%gcaAnalysisPlottingDistribPerFrame
% 
% INPUT: 
% measC: rx1 cell where r (rows) is the number of frames, each cell
% contains rx1 measurements for that frame 
% Output of  

%% Input Parser 
ip = inputParser;

ip.CaseSensitive = false;

ip.addParameter('interactive',true,@(x) islogical(x)); 
ip.addParameter('measurements',[],@(x) iscell(x)); 

% Note need to make more generic for channel wrap! 
defaultMeasDir = [MD.outputDirectory_ filesep 'MEASUREMENT_EXTRACTION' ]; 
ip.addParameter('MeasurementDirectory',defaultMeasDir,@(x) ischar(x)); 


defaultOutDir = MD.outputDirectory_; 
ip.addParameter('OutputDirectory',defaultOutDir); 



defaults{1,1} = 'filoDensityAlongVeil'; defaults{1,2} = [0,10];
defaults{2,1} = 'filoOrient'; defaults{2,2} = [0,180];
defaults{3,1} = 'filoIntensityEmbedded';defaults{3,2} = [0.5,2];
defaults{4,1} = 'filoIntensityToVeil'; defaults{4,2} = [0.5,2];
defaults{5,1} = 'filoLengthEmbedded';defaults{5,2} = [0,10];
defaults{6,1} = 'filoLengthFullActinBundle';defaults{6,2} = [0,10];
defaults{7,1} = 'filoLengthToVeil'; defaults{7,2} = [0,15];


ip.addParameter('minMaxDefaults',defaults); % defaults for me are set below 

ip.addParameter('Movie',true); 

ip.addParameter('Save',true); 

ip.addParameter('SplitMovie',true); 
ip.addParameter('SplitFrame',62); % start of the new frame 
colors(1,:) =  [  0.0314    0.1882    0.4196 ]; 
colors(2,:) =  [0.1294    0.4431    0.7098];
ip.addParameter('Color',colors); 

ip.parse(varargin{:});

%%
if ip.Results.SplitMovie 
    grpVar{1} = repmat(1,[ip.Results.SplitFrame-1,1]); 
    grpVar{2} = repmat(2,[length(ip.Results.SplitFrame:MD.nFrames_-1),1]); 
    colorGrp = vertcat(grpVar{:}); 

end 
%% go into each folder and look for a measurement 
% for now just search files - redesign so that the parameters in the
% future are more cleverly named

% might also include ylabel name and ylim for each parameter and
% read in each Descriptor directory to keep constant.

% search all descriptor parameters.
localParamFiles = searchFiles('meas_',[],[ip.Results.MeasurementDirectory filesep 'Descriptor'],1);

paramNamesC = cellfun(@(x) strrep(x,'meas_',''),localParamFiles(:,1),'uniformoutput',0);
paramNamesC = cellfun(@(x) strrep(x,'.mat',''),paramNamesC,'uniformoutput',0);

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
        selected = paramNamesC; 
       
    end    
end


 %% wrap through selections       
 for iSelect = 1:numel(selected)
    % load the measC file 
       load([localParamFiles{paramSelect(iSelect),2} filesep localParamFiles{paramSelect(iSelect),1}]);
     
       if ip.Results.Movie
           % make the directory to save the measurement movie
           outDir = [localParamFiles{paramSelect(iSelect),2} filesep selected{iSelect} '_DistributionInTime_Movie'];
           
           if ~isdir(outDir)
               mkdir(outDir) ;
           end
           nFrames = numel(measC); 
       else 
           outDir = localParamFiles{paramSelect(iSelect),2};
           nFrames = 1; 
       end
  [yLims] =  defaults{strcmpi(selected{iSelect},defaults(:,1)),2};
for iFrame = 1:nFrames 
       
       % make the boxplots
       setAxis('on');
       
     
       % convert cell to array format for boxplot plotting within matlab
       dataMat = reformatDataCell(measC);
       if ip.Results.SplitMovie 
           boxplot(dataMat,'ColorGroup',colorGrp,'outlierSize',1,'Color',ip.Results.Color,'notch','on','symbol','+'); 
       else 
           color = ip.Results.Color(:,1); 
       boxplot(dataMat,'outlierSize',1,'Color',color, 'notch','on','symbol','+' );
       end 
       hold on 
       
       frames = size(dataMat,2);
       axis([0 frames yLims(1),yLims(2)])
       
       xlabel('Time (s)');
       ylabel(selected); 
       %% MAKE More generic !!! 
       set(gca,'XTick',[100/5+1,200/5+1,300/5+1,400/5+1,500/5+1,600/5+1]); 
       set(gca,'XTickLabel',{'100','200','300','400','500','600'}); 
       
       
       if ip.Results.Movie
           yC = median(measC{iFrame});
           scatter(iFrame+0.5,yC,100,ip.Results.Color,'filled');
       end
       
       if ip.Results.Save
           if ip.Results.Movie
               saveas(gcf,[outDir filesep num2str(iFrame) '.png']);
           else
               saveas(gcf,[outDir filesep selected{iSelect} '.fig']);
               saveas(gcf,[outDir filesep selected{iSelect} '.png']);
               saveas(gcf,[outDir filesep selected{iSelect} '.eps'],'psc2');
           end
       end
      close gcf 
end % iFrame

end % iSelect


