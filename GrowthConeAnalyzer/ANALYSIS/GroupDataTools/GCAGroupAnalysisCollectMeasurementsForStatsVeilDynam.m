function [ toPlot] = GCAGroupAnalysisCollectMeasurementsForStatsVeilDynam(toPlot,varargin)
% GCAGroupAnalysisCollectMeasurementsForStatsVeilDynam

%% INPUTPARSER

%%Input check
ip = inputParser;

ip.CaseSensitive = false;
ip.addRequired('toPlot');
ip.addParameter('clearOldFields',false);
defaultOutputDirectory = pwd ;
ip.addParameter('OutputDirectory',defaultOutputDirectory);
ip.addParameter('percentileThresh',0);
ip.addParameter('splitMovie',false);
ip.addParameter('splitFrame',61);
ip.addParameter('windowChoice', 'ConstantNumber'); % Current analysis folder
% where the veil window analysis is kept ... currently have run through
% with 'ConstantNumber' or 'ProtrusionBased' with a re-initiation at frame
% 61 (5 minutes) mid movie 
ip.addParameter('umPerMin',false); % default is nm per sec 

ip.addParameter('subRoi',false)

ip.parse(toPlot,varargin{:});

%%
if ~ip.Results.splitMovie 
   % stupidly named these slightly different for the CK666... fix that here
   % 
windFold = ['protrusion_samples_' ip.Results.windowChoice '_windSize_5ReInit' num2str(ip.Results.splitFrame)];
else 
    windFold = ['protrusion_samples_' ip.Results.windowChoice '_windSize_5_ReInit' num2str(ip.Results.splitFrame)];
end 

if ip.Results.clearOldFields
    params = fieldnames(toPlot);
    params = params(~strcmpi(params,'info')); % keep the info
    toPlot = rmfield(toPlot,params);
end

params{1} = 'mednVeloc';
params{2} = 'maxVeloc';
params{3} = 'minVeloc' ;
params{4} = 'Veloc';

params{5} = 'persTime';

% params{6} = 'percentage';

analType{1} = 'protrusionAnalysis';
analType{2} = 'retractionAnalysis';
if ~ip.Results.umPerMin
    ylabel{1,1} = {'Median Velocity ' ; 'of Protrusion Event ' ; '(nm/sec)'};
    ylabel{1,2} = {'Median Velocity ' ; 'of Retraction Event ' ; '(nm/sec)'};
   
    ylabel{2,1} = {'Max Velocity ' ; 'of Protrusion Event ' ; '(nm/sec)'};
    ylabel{2,2} = {'Max Velocity ' ; 'of Retraction Event ' ; '(nm/sec)'};
    ylabel{3,1} = {'Min Velocity ' ; 'of Protrusion Event ' ; '(nm/sec)'};
    ylabel{3,2} = {'Min Velocity ' ; 'of Retraction Event ' ; '(nm/sec)'};
    ylabel{4,1} = {'Mean Velocity ' ; 'of Protrusion Event ' ; '(nm/sec)'};
    ylabel{4,2} = {'Mean Velocity ' ; 'of Retraction Event ' ; '(nm/sec)'};
    
    ylabel{5,1} = {'Persistence' ; 'of Protrusion Event ' ; '(s)'};
    ylabel{5,2} = {'Persistence' ; 'of Retraction Event ' ; '(s)'};
    
    %     ylabel{6,1} = {'Percentage of Windows' ; 'Protruding'};
    %     ylabel{6,2} = {'Percentage of Windows' ; 'Retracting'};
    %
    ylim{1} = 100; % 100 nm/sec for velocity
    ylim{2} = 100;
    ylim{3} = 100;
    ylim{4} = 100;
    ylim{5} = 60;
    ylim{6}= 1;
else
    ylabel{1,1} = {'Median Velocity ' ; 'of Protrusion Event ' ; '(um/min)'}; 
    ylabel{1,2} = {'Median Velocity ' ; 'of Retraction Event ' ; '(um/min)'};
  
    ylabel{2,1} = {'Max Velocity ' ; 'of Protrusion Event ' ; '(um/min)'};
    ylabel{2,2} = {'Max Velocity ' ; 'of Retraction Event ' ; '(um/min)'};
    ylabel{3,1} = {'Min Velocity ' ; 'of Protrusion Event ' ; '(um/min)'};
    ylabel{3,2} = {'Min Velocity ' ; 'of Retraction Event ' ; '(um/min)'};
    ylabel{4,1} = {'Mean Velocity ' ; 'of Protrusion Event ' ; '(um/min)'};
    ylabel{4,2} = {'Mean Velocity ' ; 'of Retraction Event ' ; '(um/min)'};
    
    
    ylabel{5,1} = {'Persistence' ; 'of Protrusion Event ' ; '(sec)'};
    ylabel{5,2} = {'Persistence' ; 'of Retraction Event ' ; '(sec)'};
    
    ylim{1} = 6;
    ylim{2} = 60;
    ylim{3} = 6;
    ylim{4} = 6;
    ylim{5} = 6;
    ylim{6}= 1;
end

% Collect
for iGroup = 1:numel(toPlot.info.names)
    
    projListC = toPlot.info.projList{iGroup};
    
    [nplots,~] = size(projListC);
    if ip.Results.splitMovie
        % Grouping Var1 : grouping per condition
        % create the grouping variable for pooling full group data
        % [1,(Repeated 2*nCellsProj1 times), 2(Repeated
        % 2*nCellsProj2)....[n,(Repeated 2*ncellsProjN times)]
        grpVar{iGroup} = repmat(iGroup,nplots*2,1); %
        
        % create the grouping variable for examining data per cell
        
        % Grouping Var2  :   grouping per cell
        % [1,1,2,2,3,3,...numCellsTotal,numCellsTotal]
        g = arrayfun(@(x) repmat(x,2,1),1:nplots,'uniformoutput',0);
        grpVar2{iGroup} = size(vertcat(toPlot.info.projList{1:iGroup-1}),1) + vertcat(g{:});
        
        % Grouping Var3 : grouping per treatment
        g3 = repmat([1,2],1,nplots)';
        grpVar3{iGroup} = g3 + 2*(iGroup-1);
        
    end
    
    for iAnal = 1:2
        for iParam= 1:numel(params)
            [nProjs, ~]= size(projListC);
            % initialize mat
            %dataMat = nan(7000,nProjs); % over initialize
            count = 1;
            for iProj = 1:size(projListC,1)
                folder = [projListC{iProj,1} filesep 'GrowthConeAnalyzer'];
                 
                if ip.Results.subRoi
                    toLoad1= [folder filesep 'SegmentationPackage/GCASubRegions/GC/subRegion_windows/' ip.Results.windowChoice ...
                        filesep  'EdgeVelocityQuantification_CutMovie_1' filesep 'EdgeMotion.mat'];
                    
                    toLoad2 = [folder filesep 'SegmentationPackage/GCASubRegions/GC/subRegion_windows/' ip.Results.windowChoice ...
                        filesep  'EdgeVelocityQuantification_CutMovie_2' filesep 'EdgeMotion.mat'];
                else
                    
                    toLoad1 = [folder filesep windFold filesep ...
                        'EdgeVelocityQuantification_CutMovie_1'  filesep ...
                        'EdgeMotion.mat'  ];
                    toLoad2 = [folder filesep windFold filesep ...
                        'EdgeVelocityQuantification_CutMovie_2' filesep 'EdgeMotion.mat'  ];
                end
                
                if exist(toLoad1,'file')~=0
                    first =  load(toLoad1);
                    last = load(toLoad2);
                    
                    valuesC1 =  first.analysisResults.(analType{iAnal}).total.(params{iParam}); % collect all the prot and or retract params
                    valuesC1 = valuesC1(~isnan(valuesC1));
                    
                    
                    valuesC2 = last.analysisResults.(analType{iAnal}).total.(params{iParam});
                    valuesC2 = valuesC2(~isnan(valuesC2));
                    
                    if ip.Results.percentileThresh >0
                        allVals  = [valuesC1;valuesC2];
                        thresh = prctile(allVals,ip.Results.percentileThresh);
                        valuesC1 = valuesC1(valuesC1 > thresh);
                        valuesC2= valuesC2(valuesC2 > thresh);
                        
                    end
                    
                    if ip.Results.splitMovie
                        % for now just save each as an individual cell
                        data{count} = valuesC1;
                        data{count+1} = valuesC2;
                        
                        count = count+2;
                    else
                        
                        data{count} = [valuesC1;valuesC2];
                        count = count +1;
                    end
                    
                else
                    display(['No File Found for ' projListC{iProj,1}]);
                end
            end
            
            
            dataMat = reformatDataCell(data);
            
            if ip.Results.umPerMin && iParam~=5 
                
            dataMat = dataMat.*0.06; % convert to um/min
            
            end
            
            if ip.Results.percentileThresh >0
                add = ['greaterThan' num2str(ip.Results.percentileThresh) 'th_Perc'];
            else
                add = [];
            end
            
            
            toPlot.([analType{iAnal} '_' params{iParam} '_' add]).dataMat{iGroup} = dataMat;
            toPlot.([analType{iAnal} '_' params{iParam} '_' add]).yLabel = ylabel{iParam,iAnal};
            toPlot.([analType{iAnal} '_' params{iParam} '_' add]).ylim = ylim{iParam};
            
            %
            clear data dataMat
            
        end % iParam
        
    end % iAnal
    
    
    if ip.Results.splitMovie == false
        toPlot.info.groupingPerCell = 1:size(vertcat(toPlot.info.projList{:}),1);
        toPlot.info.groupingPoolWholeMovie = toPlot.info.grouping;
    else
        toPlot.info.groupingPoolWholeMovie = vertcat(grpVar{:});
        toPlot.info.groupingPerCell= vertcat(grpVar2{:});
        toPlot.info.groupingPoolBeginEndMovie = vertcat(grpVar3{:});
    end
    
    toPlot.info.collectVeilParams= ip.Results;
    
    if ~isempty(ip.Results.OutputDirectory)
        if ip.Results.subRoi
            add = 'SubRoiGC'; 
        else 
            add = []; 
        end 
        
        if ip.Results.umPerMin; 
                add2 = 'umPerMin'; 
            else 
                add2 = []; 
        end 
            
            
        
            save([ip.Results.OutputDirectory filesep 'toPlotMeasWithVeil' add '_' add2 '.mat'],'toPlot');
        
    end
    

end % for iGroup

 
