function [ toPlot] = GCAGroupAnalysisCollectMeasurementsForStats( toPlot ,varargin)
%GCACollectData :
% INPUT: toPlot structure with field .info
%                                         .names : name of group
%                                         .projList : projList associated
%                                          with group with individual movie IDs
%                                         .grouping : vector of group
%                                          numbers
% OUTPUT: toPlot structure with all parameters that have been run added
% Note was GCACollectGroupData until 20151027
%%Input check
ip = inputParser;

ip.CaseSensitive = false;
ip.addRequired('toPlot');
ip.addParameter('clearOldFields',false);
ip.addParameter('splitMovie',false);
ip.addParameter('splitFrame', 62);  % last frame you want to include
ip.addParameter('OutputDirectory',pwd);
ip.addParameter('MeasurementFolder','SegmentationPackage/StepsToReconstructTestingGeometry20160205/GCAMeasurementExtraction_test20160221/WholeNeurite'); 
ip.addParameter('perFrame',false); % will collect the median value per frame 

ip.addParameter('filterOutlierBranchParameters',true); 
 
ip.parse(toPlot,varargin{:});
%%
nGroups = numel(toPlot.info.names);
if ip.Results.clearOldFields
    params = fieldnames(toPlot);
    params = params(~strcmpi(params,'info')); % keep the info
    toPlot = rmfield(toPlot,params);
end


for iGroup = 1:nGroups
    
    % eventually it would be nice to have a matrix detailing what was run
    % already for could just check through for the entire data set.
    
    projListC = toPlot.info.projList{iGroup}(:,1);

    nMovies = size(projListC,1);
    
    if ip.Results.splitMovie == true
        
        % Grouping Var1 : grouping per condition
        % create the grouping variable for pooling full group data
        % [1,(Repeated 2*nCellsProj1 times), 2(Repeated
        % 2*nCellsProj2)....[n,(Repeated 2*ncellsProjN times)]
        grpVar{iGroup} = repmat(iGroup,nMovies*2,1); %
        
        % Grouping Var2  :   grouping per cell
        % [1,1,2,2,3,3,...numCellsTotal,numCellsTotal]
        g = arrayfun(@(x) repmat(x,2,1),1:nMovies,'uniformoutput',0);
        grpVar2{iGroup} = size(vertcat(toPlot.info.projList{1:iGroup-1}),1) + vertcat(g{:});
        
        % Grouping Var3 : grouping per treatment
        g3 = repmat([1,2],1,nMovies)';
        grpVar3{iGroup} = g3 + 2*(iGroup-1);
        
        
        
    end % ip.Results.splitMovie
    
    % check for consistency among the parameters that were run.
    for iProj = 1:nMovies
        
        load([projListC{iProj} filesep 'GrowthConeAnalyzer'  filesep 'movieData.mat']);
        
        parameterDir = [MD.outputDirectory_ filesep ip.Results.MeasurementFolder]; 
        % for now just search files - redesign so that the parameters in the
        % future are more cleverly named
        
        % might also include ylabel name and ylim for each parameter and
        % read in each Descriptor directory to keep constant.
        
        % search all descriptor parameters.
        localParamFiles = searchFiles('meas_',[],[parameterDir filesep 'Descriptor'],1);
        
        paramNamesC = cellfun(@(x) strrep(x,'meas_',''),localParamFiles(:,1),'uniformoutput',0);
        paramNamesC = cellfun(@(x) strrep(x,'.mat',''),paramNamesC,'uniformoutput',0);
        idxOut1 = cellfun(@(x) strcmpi(x,'maxACFLagSpatial'),paramNamesC);
        idxOut2 = cellfun(@(x) strcmpi(x,'ExpressionNormalization'),paramNamesC); 
        paramNamesC(idxOut1 |idxOut2) = [];
        localParamFiles(idxOut1|idxOut2,:) = [];
        
        for iParam = 1:numel(paramNamesC)
            
            % collect the data for each parameter.
            load([localParamFiles{iParam,2} filesep localParamFiles{iParam,1}]);
            
            if ip.Results.filterOutlierBranchParameters
            
            if strcmpi(paramNamesC{iParam},'branchOrientation_2ndOrder') 
             %% Add a quick fix for some of the problems with Branch Orientation and 
                    badMeas = cellfun(@(x) x > 165,measC,'uniformoutput',0); % for now assume that any measurement larger than
                    % 165 is likely a bad orientation calculation due to
                    % the backtracing bug
                    badFrames = find(cellfun(@(x) sum(x~=0),badMeas));
                    
                    if  ~isempty(badFrames)
                        save([localParamFiles{iParam,2} filesep 'bugBranchOrientGreaterThan165.mat'],'badFrames')
                        display('badFrames found'); 
                        % filter out the problems 
                        measC =  cellfun(@(x,y) x(~y),measC,badMeas,'uniformoutput',0);
                        
                        inDir = [MD.outputDirectory_ filesep '/SegmentationPackage/StepsToReconstructTestingGeometry20160205/' ...
                            'VII_filopodiafits_geoThreshEmbed_0pt5/Channel_1'];
                        
                        
                        measToPlot{1} = 'branchOrientation_2ndOrder';
                        GCAVisualsMakeMeasurementMovieWithSub(MD,'interactive',false, ...
                            'measurements',measToPlot,...
                            'MeasurementDir',parameterDir,...
                            'InputDirectory',inDir,'frames',badFrames,'plotText',true);
                      
                         
%                          GCATroubleshootMakeMovieOfReconstructMovie(MD,'InputDirectory',...
%                                 inDir,'OutputDirectory',localParamFiles{iParam,2},'frames',badFrames);
                        
                        
                    end % badFrames 
                    
            
            end % strcmpi 'branchOrientation_2ndOrder' 
            
            if strcmpi(paramNamesC{iParam},'branchDensity_2ndOrder')
                badMeas2 = cellfun(@(x) x > 100,measC,'uniformoutput',0);
                badFrames2 = find(cellfun(@(x) sum(x~=0),badMeas2));
                if ~isempty(badFrames2)
                    save([localParamFiles{iParam,2} filesep 'bugBranchDensityGreaterThan100.mat'],'badFrames2')
                    
                    inDir = [MD.outputDirectory_ filesep '/SegmentationPackage/StepsToReconstructTestingGeometry20160205/' ...
                        'VII_filopodiafits_geoThreshEmbed_0pt5/Channel_1'];
                    
                    % filter out the problems
                    measC = cellfun(@(x,y) x(~y), measC,badMeas2,'uniformoutput',0);
                    
                    measToPlot{1} = 'branchDensity_2ndOrder';
                    
                    GCAVisualsMakeMeasurementMovieWithSub(MD,'interactive',false, ...
                        'measurements',measToPlot,...
                        'MeasurementDir',parameterDir,...
                        'InputDirectory',inDir,'frames',badFrames2,'plotText',true,...
                        'ColorByValue', false );
%                     
%                     GCATroubleshootMakeMovieOfReconstructMovie(MD,'InputDirectory',...
%                         inDir,'OutputDirectory',localParamFiles{iParam,2},'frames',badFrames2);
                    
                end
            end % branch Density 2nd order
                
            end  % filterOutlierBranchParameters
            if ip.Results.splitMovie == true
                
                % currently assumes only splitting movie in two pool these
                % values
                dataSetGroup.(paramNamesC{iParam}).valuesWholeMovie{2*(iProj-1)+1} = vertcat(measC{1:ip.Results.splitFrame});
                dataSetGroup.(paramNamesC{iParam}).valuesWholeMovie{2*(iProj-1)+2} = vertcat(measC{ip.Results.splitFrame+1:end});
                
            else
                dataSetGroup.(paramNamesC{iParam}).valuesWholeMovie{iProj} = vertcat(measC{:});
            end
            
            %% add the per frame values 
            if ip.Results.perFrame 
                if ip.Results.splitMovie == true 
                   dataSetGroup.(paramNamesC{iParam}).valuesPerFrame{2*(iProj-1)+1} = vertcat(measC{1:ip.Results.splitFrame});
                   dataSetGroup.(paramNamesC{iParam}).valuesPerFrame{2*(iProj-1)+2} = vertcat(measC{ip.Results.splitFrame+1:end});  
                else 
                    % take out empty cells 
                    %measC = measC(cellfun(@(x) ~isempty(x), measC)); 
                    idx = find(cellfun(@(x) isempty(x),measC));
                    % quick stupid fix for now
                    for i = 1:length(idx)
                        measC{idx(i)} = NaN;
                    end
                    
                    
                   
                    
                    
                    
                    %% 
                    
                    
                    valuesCell = cellfun(@(x) nanmedian(x), measC,'uniformoutput',0); % take the median per frame 
                    valuesCell = vertcat(valuesCell{:}); 
                    % for now always truncate from 1:119 
                    valuesCell = valuesCell(1:119); 
                    dataSetGroup.(paramNamesC{iParam}).valuesPerFrame{iProj} = valuesCell; 
                end 
        
            end 
            
            
        end
        
       
        
    end
    %     % collect params and reformat
    paramsAll  = fieldnames(dataSetGroup);
    %
    reformat = arrayfun(@(i) reformatDataCell(dataSetGroup.(paramsAll{i}).valuesWholeMovie),1:numel(paramsAll),...
        'uniformoutput',0);
    if ip.Results.perFrame
        reformat2 = arrayfun(@(i) reformatDataCell(dataSetGroup.(paramsAll{i}).valuesPerFrame),1:numel(paramsAll),...
            'uniformoutput',0);
        %grpPerFrame = arrayfun(@(i) 
        
        
    end
    
    for iParam = 1:numel(paramsAll)
        toPlot.(paramsAll{iParam}).dataMat{iGroup} = reformat{iParam};
        
        if  ip.Results.perFrame
            toPlot.(paramsAll{iParam}).dataMatPerFrame{iGroup} = reformat2{iParam};
            
        end
        
    end
    
    
    
    clear reformat paramsAll dataSetGroup
end  % iGroup

if ip.Results.splitMovie == false
    toPlot.info.groupingPerCell = 1:size(vertcat(toPlot.info.projList{:}),1);
    toPlot.info.groupingPoolWholeMovie = toPlot.info.grouping;
else
    toPlot.info.groupingPoolWholeMovie = vertcat(grpVar{:});
    toPlot.info.groupingPerCell= vertcat(grpVar2{:});
    toPlot.info.groupingPoolBeginEndMovie = vertcat(grpVar3{:});
    
end

%% save Data per frame 

   
  

if ~isempty(ip.Results.OutputDirectory)
    save([ip.Results.OutputDirectory filesep 'toPlotGroupMeas.mat'],'toPlot');
 
end
end