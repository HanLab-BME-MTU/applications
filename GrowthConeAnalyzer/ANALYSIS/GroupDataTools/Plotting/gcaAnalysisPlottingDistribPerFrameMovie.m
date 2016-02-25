function [filenamePNG] = gcaAnalysisPlottingDistribPerFrameMovie(MD,varargin)
%gcaAnalysisPlottingDistribPerFrame
% 
% INPUT: 
% measC: rx1 cell where r (rows) is the number of frames, each cell
% contains rx1 measurements for that frame 
% Output of  

%% Input Parser 
ip = inputParser;

ip.CaseSensitive = false;
ip.KeepUnmatched = true;
ip.addParameter('interactive',true,@(x) islogical(x)); 
ip.addParameter('measurements',[],@(x) iscell(x)); 

% Note need to make more generic for channel wrap! 

ip.addParameter('MeasurementDirectory',[],@(x) ischar(x) || isempty(x)); 


defaults{1,1} = 'filoDensityAlongVeil'; defaults{1,2} = [0,10];
defaults{2,1} = 'filoOrient'; defaults{2,2} = [0,180];
defaults{3,1} = 'filoIntensityEmbedded';defaults{3,2} = [0.5,2];
defaults{4,1} = 'filoIntensityToVeil'; defaults{4,2} = [0.5,2];
defaults{5,1} = 'filoLengthEmbedded';defaults{5,2} = [0,10];
defaults{6,1} = 'filoLengthFullActinBundle';defaults{6,2} = [0,10];
defaults{7,1} = 'filoLengthToVeil'; defaults{7,2} = [0,15];
defaults{8,1} = 'branchLength_2ndOrder'; defaults{8,2} = [0,15]; 
defaults{9,1} = 'branchOrientation_2ndOrder'; defaults{9,2} = [0 180];
defaults{10,1} = 'veilStemThickness'; defaults{10,2} = [0,10]; 

ip.addParameter('minMaxDefaults',defaults); % defaults for me are set below 

ip.addParameter('Movie',false); 

ip.addParameter('Save',true); 

ip.addParameter('SplitMovie',false); 
ip.addParameter('SplitFrame',62); % start of the new frame 

ip.addParameter('Color',[]); 
% 
ip.addParameter('visible','off'); 
ip.addParameter('CollectDirectory',[]); % 

ip.addParameter('ColorByKS',true); 

ip.parse(varargin{:});

%%
if isempty(ip.Results.MeasurementDirectory)
    measDir = [MD.outputDirectory_ filesep 'MEASUREMENT_EXTRACTION' ];
else
    measDir = ip.Results.MeasurementDirectory;
end

if isempty(ip.Results.CollectDirectory)
    collectDir = [measDir filesep ...
        'Descriptor' filesep 'Filopodia' filesep 'MeasurementsInTime'];
else
    collectDir = ip.Results.CollectDirectory;
end

if isempty(ip.Results.Color)
    if ip.Results.SplitMovie
        colors(1,:) =  [  0.0314    0.1882    0.4196 ];
        colors(2,:) =  [0.1294    0.4431    0.7098];
    else
        colors(1,:) = [0 0 0]; % black
        colors(2,:) = [0 0 1]; % blue
        colors(3,:) = [1 0 0]; % red
        
    end
else
    colors = ip.Results.colors;
end

if ip.Results.SplitMovie 
    grpVar{1} = repmat(1,[ip.Results.SplitFrame-1,1]); 
    grpVar{2} = repmat(2,[length(ip.Results.SplitFrame:MD.nFrames_-1),1]); 
    colorGrp = vertcat(grpVar{:}); 
end 



% movieID 
ID = helperGCACreateID(MD); 

% make the collection directory
if ~isdir(collectDir);
    mkdir(collectDir);
end

if ip.Results.ColorByKS
    forName = 'ColoredByKS'; 
else 
    forName = ''; 
end 
    

%% go into each folder and look for a measurement 
% for now just search files - redesign so that the parameters in the
% future are more cleverly named

% might also include ylabel name and ylim for each parameter and
% read in each Descriptor directory to keep constant.

% search all descriptor parameters.
localParamFiles = searchFiles('meas_',[],[measDir filesep 'Descriptor'],1);



paramNamesC = cellfun(@(x) strrep(x,'meas_',''),localParamFiles(:,1),'uniformoutput',0);
paramNamesC = cellfun(@(x) strrep(x,'.mat',''),paramNamesC,'uniformoutput',0);

% For now take out spatial autocorrelation 
exclude = cellfun(@(x) strcmpi(x,'maxACFLagSpatial'),paramNamesC);
paramNamesC = paramNamesC(~exclude); 
localParamFiles = localParamFiles(~exclude,:); 

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
        paramSelect = 1:numel(paramNamesC); 
    end    
end


 %% wrap through selections       
 for iSelect = 1:numel(selected)
     
     % load the measC file 
       load([localParamFiles{paramSelect(iSelect),2} filesep localParamFiles{paramSelect(iSelect),1}]);
       statFlag = true; 
       % test for an empty cell at the end and remove if necessary
       toRemove = cellfun(@(x) isempty(x), measC);
       if sum(toRemove) >0
           if find(toRemove) == numel(measC)
               % just remove the final values
               measC(toRemove) = [];  
           else
             statFlag = false ;
               
           end
       end 
      % Reinitiate Colors
      colorsC = colors; 
       if ip.Results.ColorByKS
           
           if statFlag
               
               ksResults = arrayfun(@(i) kstest2(measC{1},measC{i}),1:numel(measC)) ;
               direction = arrayfun(@(i) sign(nanmedian(measC{i})-nanmedian(measC{1})),1:numel(measC)) ;
               colorGrp = direction.*ksResults;
               colorGrp(colorGrp==1) = 3; 
               colorGrp(colorGrp==0) = 1; 
               colorGrp(colorGrp==-1) = 2; 
            
               
               % modify the colorMap based on colors 
             
               idx = unique(colorGrp);
               colorsC = colorsC(idx,:) ; 
               % modify based on colorgroups 
               
              
           else
               display(['Cannot Calculate Stats for ' ...
                    localParamFiles{paramSelect(iSelect),2} filesep localParamFiles{paramSelect(iSelect),1}]);
           end
       end
           
       
       
     
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
       
       idxD = strcmpi(selected{iSelect},defaults(:,1)); 
       if sum(idxD) ~=0
           
           
           [yLims] =  defaults{strcmpi(selected{iSelect},defaults(:,1)),2};
       else
           yLims = [];
           
       end
for iFrame = 1:nFrames 
       
       % make the boxplots
       setAxis(ip.Results.visible,0.5,14);
          
       % convert cell to array format for boxplot plotting within matlab
       dataMat = reformatDataCell(measC);
       
       if (size(dataMat,2) >1 && size(dataMat,1) >1) % check if value per frame a distribution
           if ip.Results.SplitMovie || ip.Results.ColorByKS
               if size(dataMat,2) > length(colorGrp)
                   dataMat = dataMat(:,1:length(colorGrp)); 
               end 
               if length(colorGrp) > size(dataMat,2)
                   colorGrp = colorGrp(1:size(dataMat,2)); 
               end 
               boxplot(dataMat,'ColorGroup',colorGrp,'outlierSize',1,'Color',colorsC,'notch','on','symbol','+');
           else
               colorsC = colors(1,:);
               boxplot(dataMat,'outlierSize',1,'Color',colorsC, 'notch','on','symbol','+' );
           end
       else
           color1 = ip.Results.Color(1,:);
           color2 = ip.Results.Color(2,:);
           if ip.Results.SplitMovie
               
               scatter(1:length(dataMat(1:ip.Results.SplitFrame-1)),dataMat(1:ip.Results.SplitFrame-1),50,color1,'filled');
               hold on
               scatter(ip.Results.SplitFrame:length(dataMat),dataMat(ip.Results.SplitFrame:length(dataMat)),50,color2,'filled');
           else
               
               scatter(1:length(dataMat),dataMat,50,color1,'filled');
           end
           
       end
       
          
       hold on 
       
       frames = size(dataMat,2);
       if ~isempty(yLims)
       axis([0 frames yLims(1),yLims(2)])
       end 
       xlabel('Time (s)');
       ylabel(selected{iSelect}); 
       
        if ip.Results.SplitMovie 
            if ~(isempty(yLims))
            ymin = yLims(1); 
            ymax = yLims(2);
            else 
                ymin = min(dataMat(:)); 
                ymax = max(dataMat(:)); 
            end 
            x = ip.Results.SplitFrame-0.5; 
           line([x ,x],[ymin,ymax],'linewidth',2,'color',colors(1,:)); 
       end 
       
       
       %% MAKE More generic !!! 
       set(gca,'XTick',[100/5+1,200/5+1,300/5+1,400/5+1,500/5+1,600/5+1]); 
       set(gca,'XTickLabel',{'100','200','300','400','500','600'}); 
       
      
       
       if ip.Results.Movie
           yC = median(measC{iFrame});
           scatter(iFrame+0.5,yC,100,ip.Results.Color(1,:),'filled');
       end
       
       if ip.Results.Save
           if ip.Results.Movie
               saveas(gcf,[outDir filesep num2str(iFrame) '.png']);
           else
               
               filenamePNG = [' ' outDir filesep selected{iSelect} '.png' ' ']; % for montaging 
               saveas(gcf,[outDir filesep ID '_' selected{iSelect} forName '.fig']);
               saveas(gcf,[outDir filesep ID '_' selected{iSelect} forName '.png']);
               saveas(gcf,[outDir filesep ID '_' selected{iSelect} forName '.eps'],'psc2');
               
              
                   
                   
                   saveas(gcf,[collectDir filesep ID '_' selected{iSelect} forName '.fig']); 
                   saveas(gcf,[collectDir filesep ID '_' selected{iSelect} forName '.png']);
                   saveas(gcf,[collectDir filesep ID '_' selected{iSelect} forName '.eps'],'psc2');
                   
            
               
           end
       end
      close gcf 
end % iFrame

end % iSelect


