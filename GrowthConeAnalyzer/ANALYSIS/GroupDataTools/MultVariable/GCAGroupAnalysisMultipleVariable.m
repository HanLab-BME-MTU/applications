function [ output_args ] = GCAGroupAnalysisMultipleVariable(toPlot,varargin)
% GCAGroupAnalysisMultipleVariable : A working function that compiles different
% options for exploring the multivariable neurite data (predictors) and its relationship
% to a response variable such as elongation.
% Current Options:
% MDS (multi-D scaling)
% PCA
% Stepwise Mult Var Regression:
%% REQUIRED
%
%  toPlot: (REQUIRED)  group structure for groupAnalysis with the collected data
%                      for each field
%
%% PARAMETERS
%
%    'interactive' (PARAM) : logical
%       if true the user can pick which measurements to include in the
%       analysis
%
%
%% Check Input
ip = inputParser;
ip.CaseSensitive = false;
% Check input
ip.addRequired('toPlot');

% Feature Selection Options
ip.addParameter('interactive',true);

% Regression Analyses (Will hopefully eventually use for feature selection)
ip.addParameter('SWMultReg',false); % flag to perform multiple variable regression
ip.addParameter('lasso',false);

% Function for movie stat
ip.addParameter('perNeuriteStatistic','nanmedian'); % default is to take the
% median value of the measurement calculated over the entire movie.

% InputOutput
defaultOut = pwd;
ip.addParameter('OutputDirectory',defaultOut);

% Plots Types
ip.addParameter('plotByGroupElongation',false); % will make 3-D scatter plots by group
ip.addParameter('perFrame',false); 
ip.addParameter('plotByGroupColorID',false); 
ip.addParameter('plotByGroupPic',false); 


ip.addParameter('plotCentroids',true); 


% D Reduction Visualizations
ip.addParameter('PCA',false);
ip.addParameter('MDS',true);
ip.addParameter('Run', true); % will search for a 
%ip.addParameter('CriterionMDS','stress'); % default is nonmetric

ip.addParameter('cluster',true);
ip.addParameter('clustNumResponse',3);
ip.addParameter('testResponseClust',false); 

ip.addParameter('DistMetrics',false); 

ip.addParameter('ColorTrajByTime',true); 
ip.addParameter('MakeMovie',false); 

ip.parse(toPlot,varargin{:});
%% Set up
varNames  = fieldnames(toPlot);
varNames = varNames(~strcmpi(varNames,'info')) ;

projList = vertcat(toPlot.info.projList{:});
if size(projList,2)==2
    
    obsNames = projList(:,2);
else
    obsNames = cellfun(@(x) helperGCACreateID(x),projList(:,1),'uniformoutput',0);
    
end

% ip perframe for each of the observation names repeat by 120 
if ip.Results.perFrame
    % add number to obsNames per frame
    for iMovie = 1:size(projList,1);
        obsNamePerFrame = arrayfun(@(x) [obsNames{iMovie} '_' num2str(x,'%02d')],1:119,'uniformoutput',0);
        obsNamePerFrameAll{iMovie} = obsNamePerFrame'; 
    end
    obsNamePerFrameAll = obsNamePerFrameAll'; 
    obsNames = vertcat(obsNamePerFrameAll{:});
end

% if ip.Results.perFrame
%     for i = 1:120 
%     obsNamesGrp = arrayfun(@(x) [obsNames{x} num2str(x,'%02')]); 
%     end
%     
% end
nGroups = numel(toPlot.info.names);
grouping = toPlot.info.grouping;

if ip.Results.Run 


%% User can define which measurements(variables) to include if interactive turned on


if ip.Results.interactive
    paramSelect  = listSelectGUI(varNames,[],'move');
    varNames = varNames(paramSelect);
end
nVars = numel(varNames);
varNamesString = strrep(varNames, '_' ,'');

%% plot by frame : note this is going to need to be a different structure
% for now just add a quick and dirty grouping variable for per Frame 
% make a single vector labeling each trajectory for each group (just use
% the first varName for now 
if ip.Results.perFrame
    % check if the size of all the parameters is the same (same number of
    % frames)
    for iGroup = 1:nGroups
        
        framesPerParam = arrayfun(@(x) size(toPlot.(varNames{x}).dataMatPerFrame{iGroup},1),1:nVars);
        nFrames = length(unique(framesPerParam));
        if nFrames~=1
            msg = ['Number of Frames not Consistent for All Variables for Group ' toPlot.info.names{iGroup}];
            error(msg);
        else
            nMovies = size(toPlot.info.projList{iGroup},1);
            grpingPerFrame = arrayfun(@(x) repmat(x,framesPerParam(1),1),1:nMovies,'uniformoutput',0);
         
            
            
            grpingTrajIDsAllGroups{iGroup} = vertcat(grpingPerFrame{:});
            
            
        end
    end
    grpIDs = arrayfun(@(x) repmat(x,framesPerParam(1).*size(toPlot.info.projList{x},1),1),1:nGroups,'uniformoutput',0);
    grpingGroupIDsAllGroups = vertcat(grpIDs{:}); 
    grpingTrajIDsAllGroups = vertcat(grpingTrajIDsAllGroups{:}); 
end









%% Collect and normalize the measurement data : Default is to take a median value for each cell

dataDir = [ip.Results.OutputDirectory filesep 'DataFiles'];
if ~isdir(dataDir)
    mkdir(dataDir);
end

save([dataDir filesep 'varNames'],'varNames'); 

for iGroup  = 1:numel(toPlot.info.names)
    
    
    if ip.Results.perFrame
        % make sure all measurements have the same number of frames ...
        
        % cellfun(@(x) size(
        
        % combine all neurites for each measurement, data are now a rx1
        dataInt = arrayfun(@(x) toPlot.(varNames{x}).dataMatPerFrame{iGroup}(:),1:nVars,'uniformoutput',0);
        
        dataMatC = horzcat(dataInt{:});
        
    else
        
        % transform data so that each row is a neurite and each column is a
        % set of measurements- this will be a per cell variable.
        perNeuriteStat = str2func(ip.Results.perNeuriteStatistic);
        dataMatC = arrayfun(@(x) perNeuriteStat(real(toPlot.(varNames{x}).dataMat{iGroup}),1)',1:nVars,'uniformoutput',0);
    
        % put all the measurement values together (median measurement value for each neurite per column)
        % r indicates the neurite number;
        dataMatC = horzcat(dataMatC{:});
    end 
    % dataMatAllMeas is a rxc matrix where r is the number of neurites tested
    % and c is the number of measurements
    dataMatAllMeas{iGroup} = dataMatC ; % this will keep them in a cell by KD condition
    
end

% try once by normalizing each measurement distribution
% using the pooled population regardless of condition (ie control + KD)
% this would be mainly to make sure the measurements are on the same scale.
dataMatAllGroupsMeas = vertcat(dataMatAllMeas{:});

%dataFinal = dataMatAllGroupsMeas;
dataFinal = zscoreForNaNs(dataMatAllGroupsMeas); % this should by default normalize along the column

% put back into the full group into a cell 
dataFinalAll = arrayfun(@(x) dataFinal(grouping == x,:),1:nGroups,'uniformoutput',0); 


if ip.Results.plotCentroids
    centroidsOfGroups = arrayfun(@(x) nanmean(dataFinal(grouping==x,:),1),1:nGroups,'uniformoutput',0); % find the centroid of the features 
    centroidsOfGroups = vertcat(centroidsOfGroups{:}); 
    dataFinal = [dataFinal;centroidsOfGroups];
    % make the centroid points grouping == to zero 
    groupingCent = zeros(nGroups,1); 
    grouping = [grouping;groupingCent]; 
    obsNamesCent = cellfun(@(x) ['Centroid ' x] , toPlot.info.names,'uniformoutput',0); 
    obsNamesZ = [obsNames;obsNamesCent']; 
   
else 
    obsNamesZ = obsNames; 
end 



% Save the raw and normalized collected data in a csv.
% NOTE: need to change to table eventually to make compatible with later
% matlab versions!

% if ip.Results.perFrame
%     dataPreName = mat2dataset(dataMatAllGroupMeas,'varNames',varNames); 
% else 
dataPreNorm = mat2dataset(dataMatAllGroupsMeas,'ObsNames',obsNames,'varNames',varNames);
% end 
save([dataDir filesep ['OriginalValues' ip.Results.perNeuriteStatistic 'PerNeuriteMovie']],'dataPreNorm');

export(dataPreNorm,'file',[dataDir filesep 'OriginalValuesMedianPerNeuriteMovie.csv']);
%writetable(tDataRaw,[ip.Results.OutputDirectory filesep 'OriginalValuesMedianPerNeuriteMovie.csv'],dataOriginal);

dataZScores = mat2dataset(dataFinal,'ObsNames',obsNamesZ,'varNames',varNames);
export(dataZScores,'file',[dataDir filesep 'dataZScores.csv']);
save([dataDir filesep 'dataZScores'],'dataZScores');
else % load the data 
     
    
    
end 
%% Collect the Outgrowth Data, Create Mapper, and Save Colorbar

outgrowth = gcaCollectOutgrowthDeltasPerGroup(toPlot);
outgrowth = vertcat(outgrowth{:});
outgrowth = outgrowth./10;

 
    fsFigure(0.75,'visible','off');
    test = -2:2;
    imagesc(test);
    
    cmap = brewermap(128,'RdBu');
    cmap = flip(cmap,1);
    % get the average velocity values for all the windows in the current frame
    % and assign a color based on the the mapper.
    plotValues = outgrowth;
    
    mapper=linspace(-2,2,128)'; % lower values are red in the brewermap
    D=createDistanceMatrix(plotValues,mapper);
    [sD,idxCMap]=sort(abs(D),2);
    colormap(cmap);
    colorbar
    saveas(gcf,[dataDir filesep 'colorbar.fig']);
    saveas(gcf,[dataDir filesep 'colorbar.eps'],'psc2');

    close gcf



if ip.Results.cluster
    clusterDir = [ip.Results.OutputDirectory filesep 'Cluster'];
    if ~isdir(clusterDir)
        mkdir(clusterDir);
    end
    
    %% test response clusters
    if ip.Results.testResponseClust
        
        for i = 1:5
       
        myfunc = @(X,K)(kmeans(X, K, 'replicates',50)); 
            
            
        eval = evalclusters( plotValues, myfunc,'DaviesBouldin','klist',1:10);
        evalControl = evalclusters(plotValues(grouping ==1),myfunc,'DaviesBouldin','klist',1:10);
        
        setAxis('on')
        plot(1:10,evalControl.CriterionValues,'color','k');
        
        hold on
        plot(1:10,eval.CriterionValues,'color','r');
        xlabel('Cluster Number');
        ylabel('DaviesBouldin Index');
        nControl = length(plotValues(grouping==1));
        nCells = length(plotValues);
        legend({['Control N = ' num2str(nControl) ], ['All Conditions  N = ' num2str( nCells  )] }  );
        legend('boxoff');
        scatter(1:10,evalControl.CriterionValues,50,'k','filled');
        scatter(1:10,eval.CriterionValues,50,'r','filled');
        saveas(gcf,[clusterDir filesep 'DBIndices' num2str(i) '.fig']);
        saveas(gcf,[clusterDir filesep 'DBIndices' num2str(i) '.png']);
        close gcf
        end 
       
    end % if testResponseClust
   %scatter(eval
%     colorClust{1} = [0 0 1]; % low blue
%     colorClust{2} = [1 0 0]; % high red
%     colorClust{3} = [0 1 0]; % mid
%      
     %cMap =  brewermap(ip.Results.ClustNumResponse,'set1'); 
     cMap(:,1) = [0,0,1]; 
     cMap(:,2) = [0,1,0]; 
     cMap(:,3) = [1,0,0]; 
    % Perform k-means clustering (All Values)
    [idx,cCenters]  = kmeans(plotValues,ip.Results.clustNumResponse,'replicates',50);
    %     indexMin = find(cCenters == min(cCenters));
    %     indexMax = find(cCenters == max(cCenters));
    [values,indices] = sort(cCenters);
   % indices
   sortedIdx = zeros(length(idx),1);
   for i = 1:length(indices)
       
       sortedIdx(idx==indices(i)) =i;
   end
   groupingCluster = sortedIdx;
   if ip.Results.plotCentroids
       pad = zeros(nGroups,1); 
       groupingCluster = [groupingCluster;pad]; 
       
   end 
    
%     %  % make it such that 2 is always high and 1 is always the low cluster
%     sortedIdx = zeros(length(idx),1);
%     sortedIdx(idx==indexMin) = 1;
%     sortedIdx(idx==indexMax) = 2;
%     groupingCluster = sortedIdx;
    


%% add the plot values and the cluster to save 
end

%% Perform Stepwise linear regression
if ip.Results.SWMultReg
    SWDir = [ip.Results.OutputDirectory filesep 'StepWiseMultReg'];
    if ~isdir(SWDir)
        mkdir(SWDir);
    end
    [coef,se,pval,inmodel,stats,nextstep,history] = stepwisefit(dataMatAllGroupsMeas,outgrowth);
    % there is also stepwiselm creates an object!
    resultsStepWiseFitAll = [num2cell(coef(inmodel)) varNamesString(inmodel)];
    
    save([SWDir filesep 'StepWiseMultRegResultsAllData.mat'],'coef','se','inmodel','stats','nextstep','history','resultsStepWiseFitAll');
    
    
    dataMatControl = dataMatAllGroupsMeas(1:size(toPlot.info.projList{1}),:);
    outgrowthControl = outgrowth(1:size(toPlot.info.projList{1}));
    
    [coefCon,seCon,pvalCon,inmodelCon,statsCon,nextstepCon,historyCon] = stepwisefit(dataMatControl,outgrowthControl);
    
    resultsStepWiseFitCon = [num2cell(coefCon(inmodelCon)) varNamesString(inmodelCon)];
    
    
    save([SWDir filesep 'StepWiseMultRegResultsAllControl.mat'],'coefCon','seCon','inmodelCon','statsCon','nextstepCon','historyCon','resultsStepWiseFitCon');
    
end
%% Perform Lasso
if ip.Results.lasso
    lassoDir= [ip.Results.OutputDirectory filesep 'lasso'];
    if ~isdir(lassoDir)
        mkdir(lassoDir);
    end
    
    [coef,fitInfo] = lasso(dataMatAllGroupsMeas,outgrowth);
    %     %% there is also stepwiselm creates an object
    %
    %
    %     dataMatControl = dataMatAllGroupsMeas(1:size(toPlot.info.projList{1}),:);
    %     outgrowthControl = outgrowth(1:size(toPlot.info.projList{1}));
    %
    %     [coefCon,seCon,pvalCon,inmodelCon,statsCon,nextstepCon,historyCon] = stepwisefit(dataMatControl,outgrowthControl);
    %
    %     resultsStepWiseFitCon = [num2cell(coefCon(inmodelCon)) varNamesString(inmodelCon)];
    %
    %
    %        save([SWDir filesep 'StepWiseMultRegResultsAllControl.mat'],'coefCon','seCon','inmodelCon','statsCon','nextstepCon','historyCon','resultsStepWiseFitCon');
    %
    %
    %
    % %     [coefCon,seCon,pvalCon,inmodelCon,statsCon,nextStepCon,historyCon] = stepwisefit(
    % %     save([SWDir filespe 'StepWiseMultRegResultsControlOnly.mat']
    
end

%% Plot tSNE: test

% tsneVals = tsne(dataFinal,[],2,0.05);
%
% colors = vertcat(toPlot.info.color{:});
% grouping =toPlot.info.grouping;
% setAxis('on')
% hold on
% arrayfun(@(x) scatter(tsneVals(grouping ==x,1),tsneVals(grouping == x,2),50,colors(x,:),'filled'),1:size(colors,1));
% saveas(gcf,'tSNEPlot.fig');
% close gcf

%% Norm II :NOTE Do NOT use this normalization this will normalize each population per condition which is not
%  really what we want.
% Normalize each measurement distribution to zero mean and standard deviation of
% 1 for each condition group- these are when each are normalized
% independently - likely not as appropriate
% dataMatZ = cellfun(@(x) zscoreForNaNs(x),dataMatAllMeas,'uniformoutput',0);
% dataFinal = vertcat(dataMatZ{:});


%% Perform the MDS multi dimensional scaling
if ip.Results.MDS
    % make MDS file
    outMDS = [ip.Results.OutputDirectory filesep 'MDS'];
    if ~isdir(outMDS);
        mkdir(outMDS);
    end
    
    
    % rows are observations, columns variables
    % compute pairwise distances
    dissimilarities = pdist(dataFinal);
    z = squareform(dissimilarities);
    
    
    if ip.Results.plotCentroids
        % find the closest point to each centroid
        % get all the rows that correspond to the centroids and all the
        % columns that correspond to each group
        nCells = size(vertcat(toPlot.info.projList{:}),1); 
        idxClosest  =   arrayfun(@(x)  find(z(nCells+x,grouping==x)==min(z(nCells+x,grouping==x))),1:nGroups,'uniformoutput',0); % rows will be
       
        idxFurthest = arrayfun(@(x) find(z(nCells+x,grouping==x)== max(z(nCells+x,grouping ==x))),1:nGroups,'uniformoutput',0); 
        IDsFurthest = arrayfun(@(x) toPlot.info.projList{x}(idxFurthest{x},:),1:nGroups,'uniformoutput',0);  
        %arrayfun(@(x) find(distMatCents{x,:} == max(i,:)),
        IDsFurthest = vertcat(IDsFurthest{:}); 
        
        [values,idxSort] = arrayfun(@(x) sort(z(nCells+x,grouping==x)),1:nGroups,'uniformoutput',0); 
        IDsClosest = arrayfun(@(x) toPlot.info.projList{x}(idxClosest{x},:),1:nGroups,'uniformoutput',0);  
        %arrayfun(@(x) find(distMatCents{x,:} == max(i,:)),
        IDsClosest = vertcat(IDsClosest{:}); 
        
    end
    
    
        test = isnan(z(:,1)); 
        
        if ip.Results.perFrame
            if sum(test)~=0
                display('NaN values found: Check Measurements'); 
                z(test,:) =[];
                z(:,test) = [];
                grpingGroupIDsAllGroups(test) =[]; 
                grpingTrajIDsAllGroups(test) = []; 
                removed = obsNames(test); 
                %removed = [removed num2cell(dataFinal(test,:))]; 
                valuesRemoved = dataFinal(test,:); 
                dataRemoved = mat2dataset(valuesRemoved,'ObsNames',removed,'varNames',varNames);
                export(dataRemoved,'file',[dataDir filesep 'dataRemove_NaN.csv']);
                save([dataDir filesep 'dataRemoved'],'dataRemoved');
                       
                % Take out the rest  
                obsNames(test) = [];
                dataFinal(test,:) = []; 
                
                % save the new data structure 
                dataClean = mat2dataset(dataFinal,'ObsNames',obsNames,'varNames',varNames); 
                export(dataClean,'file',[dataDir filesep 'zScoredDataClean_NoNaN.csv']); 
                
                % redo the dissimilarity matrix 
                dissimilarities = pdist(dataFinal);
                squareFormDis = squareform(dissimilarities); % save the dissimilarities in square form 
                save([dataDir filesep 'squareFromDis.mat'],'squareFormDis'); 
                           
            end
            
        else
            
            test = isnan(z(:,1));
            if sum(test)~=0
                z(test,:) =[];
                z(:,test) = [];
                
                
                
                grouping(test) = [];
                projListAll = vertcat(toPlot.info.projList{:});
                removed = projListAll(test,2);
                cellfun(@(x) display(['Removing cell' x 'due to NaN: Check Measurements']),removed);
                
                
                
                projListAll(test,:) = [];
                idxCMap(test,:) = [];
                plotValues(test) =[];
                
                
            end
        end  
        
     
        
        
        
        
      
    %      z= squareform(dissimilarities);
    
    
    % Test in case NaN there doesn't appear to be any NaN any more = should
    % treat as empty anyway.
    
    % test = isnan(z(:,1));
    % if sum(test)~=0
    %     z(test,:) =[];
    %     z(:,test) = [];
    %     grouping(test) = [];
    %     projListAll = vertcat(toPlot.info.projList{:});
    %     removed = projListAll(test,2);
    %     cellfun(@(x) display(['Removing cell' x 'due to NaN: Check Measurements']),removed);
    %     projListAll(test,:) = [];
    %     idxCMap(test,:) = [];
    %     plotValues(test) =[];
    %
    %
    % end  % test
    % for now do two criterion
    criterion{1} = 'stress';
%     criterion{2} = 'metricstress';
%     criterion{3} = 'sammon';
%     criterion{4} = 'strain';
    
    for iCrit = 1:numel(criterion)
        
        cOutMDS =  [outMDS filesep criterion{iCrit} ];
        if ~isdir(cOutMDS)
            mkdir(cOutMDS)
        end
        
        cOutMDSDia = [cOutMDS filesep 'Diagnostic'];
        
        if ~isdir(cOutMDSDia)
            mkdir(cOutMDSDia);
        end
        
        % check for .mat file
        if exist([cOutMDSDia filesep 'MDSResults' criterion{iCrit} '.mat'],'file')==0;
            
            display('Running MDS'); 
            % perform the nonclassical multi-dimensional scaling
            [ y, stress, disparities] = mdscale(dissimilarities,2,'Criterion',criterion{iCrit});
            % if strcmpi(criterion{iCrit},'stress');
            save([cOutMDSDia filesep 'MDSResults' criterion{iCrit}],'y','stress','disparities','dissimilarities');
            display('Finished MDS'); 
            %else
            %   save([cOutMDSDia filesep 'MDSResults' criterion{iCrit}],'y','stress');
            %end
        else
            display(['File for MDS ' criterion{iCrit} ' Found: Loading']);
            load([cOutMDSDia filesep 'MDSResults' criterion{iCrit} '.mat']);
        end
        
           
        if ip.Results.perFrame
            axisLims = [prctile(y(:,1),0.5),prctile(y(:,1),99.5),prctile(y(:,2),0.5),prctile(y(:,2),99.5)];
        else
            axisLims = [floor(min(y(:,1))-1),ceil(max(y(:,1)))+1,floor(min(y(:,2))-1),ceil(max(y(:,2))+1)];
        end
        %% Plot shepards plot
%         setAxis('on');
%         distances = pdist(y);
%         [dum,ord] = sortrows([disparities(:) dissimilarities(:)]);
%         plot(dissimilarities,distances,'bo', ...
%             dissimilarities(ord),disparities(ord),'r.-', ...
%             [0 25],[0 25],'k-')
%         %transformed values of the original dissimilarities
%         
%         xlabel('Dissimilarities')
%         ylabel('Distances/Disparities')
%         legend({'Distances' 'Disparities' '1:1 Line'},...
%             'Location','NorthWest');
%         saveas(gcf,[cOutMDSDia filesep 'shepardsPlot.fig']);
%         saveas(gcf,[cOutMDSDia filesep 'shepardsPlot.png']);
%         
%         close gcf
        
        
        %% Plot in MDS space By Group
        
        % set up figure
        setAxis('on')
        if ip.Results.perFrame
            
            % make a folder
            grpDir = [ip.Results.OutputDirectory filesep 'MDSTrajectoryPlotsPerGroup'];
            if ~isdir(grpDir);
                mkdir(grpDir);
            end
            
            for iGroup = nGroups-1:nGroups-1
                trajDir = [grpDir filesep toPlot.info.names{iGroup}];
                if ~isdir(trajDir)
                    mkdir(trajDir);
                end
                                   
                nTrajs = size(toPlot.info.projList{iGroup},1);
                cmapCurrent = repmat(toPlot.info.colorShades{iGroup},3,1);
                
                projListC = toPlot.info.projList{iGroup};
                
                for iTraj = 1:nTrajs
                setAxis('off')
                % plot all results in gray
                scatter(y(:,1),y(:,2),50,[.5,.5,.5],'filled');
                hold on
                
                currentTraj = find(grpingTrajIDsAllGroups==iTraj & grpingGroupIDsAllGroups ==iGroup);
                if ~isempty(currentTraj)
                %if ~isempty(currentTraj)
%                 imgStart = imread([projListC{iTraj,1} filesep 'GrowthConeAnalyzer' filesep 'VisualizationOverlays' filesep 'Raw'...
%                     filesep '001.png']);
%                 imgEnd = imread([projListC{iTraj,1} filesep 'GrowthConeAnalyzer' filesep 'VisualizationOverlays' filesep 'Raw'...
%                     filesep '118.png']);
%                 image([y(currentTraj(1),1)-0.01, y(currentTraj(1),1)+0.01],[y(currentTraj(1),2)-0.01,y(currentTraj(1),2)+0.01],imgStart);
%                 image([y(currentTraj(end),1)-0.01, y(iTraj,1)+0.01],[y(iTraj,2)-0.01,y(iTraj,2)+0.01],imgEnd);
%                 
                
                
               
                %for iTraj = 1:nTrajs
                    
                    if ip.Results.ColorTrajByTime
                        % make individual file for each 
                         
                        %cmapTime = jet(framesPerParam(1));
                        cmapTimeBeforeDrug  = repmat(toPlot.info.colorShades{10}(end,:),[61,1]) ; ... 
                            
                        cmapTimeAfterDrug = repmat(toPlot.info.colorShades{10}(5,:),[61,1]); 
                        
                        cmapTime= [cmapTimeBeforeDrug ; cmapTimeAfterDrug];
                        
                        
                        %if ip.Results.ColorTrajByTimeMovie 
                            
                            name = obsNames{currentTraj(1)}; 
                            name = strrep(name,'_','/');
                            name = upDirectory(name,1); 
                            name = strrep(name,'/',' '); 

                            title(['Trajectory ' name], 'FontSize',14); 
                            nameUnder =  strrep(name,' ' ,'_');
                            
                            
                            
                            perTrajDir = [trajDir filesep 'Traj_' nameUnder ];
                            if ~isdir(perTrajDir)
                                mkdir(perTrajDir)
                            end
                            
                            
                            if ip.Results.MakeMovie
                                frames = length(currentTraj);
                                for i = 1:length(currentTraj)
                                    
                                    arrayfun(@(i) plot(y(currentTraj(i):currentTraj(i+1),1),y(currentTraj(i):currentTraj(i+1),2),...
                                        'color',cmapTime(i,:)),1:length(currentTraj)-1);
                                    scatter(y(currentTraj(i),1),y(currentTraj(i),2),50,cmapTime(i,:),'filled')
                                    axis(axisLims);
                                    
                                    xlabel('MDS1');
                                    ylabel('MDS2');
                                    
                                    
                                    
                                    saveas(gcf,[perTrajDir filesep num2str(i,'%03d') '.png']);
                                    
                                    
                                    close gcf
                                    setAxis('off')
                                    % plot all results in gray
                                    scatter(y(:,1),y(:,2),50,[.5,.5,.5],'filled');
                                    title(['Trajectory ' name 'Time ' num2str(i*5) 'sec'], 'FontSize',14);
                                    hold on
                                    
                                end
                                
                                % setUpTheMontage
                                inputFolders{1} = perTrajDir; 
                                inputFolders{2} = [projListC{iTraj} filesep ... 
                                    'GrowthConeAnalyzer/SegmentationPackage/StepsToReconstructTestingGeometry20160205/' ...
                                    'GCAVisuals/Descriptor/Filopodia/Validation/filoLength/ForMainMovie_Measurement_Movie']; 
                                
                                outputFolder = [perTrajDir filesep 'Montages']; 
                                if ~isdir(outputFolder)
                                    mkdir(outputFolder);
                                end 
                                
                              GCAVisualsMontagingMovieNoMD(frames,inputFolders,outputFolder); 
                            end
                            
                            
                            
                            
                            
                            
                            % Set up Montage
                            imageFolder  = ['GrowthConeAnalyzer/SegmentationPackage/' ...
                                'StepsToReconstructTestingGeometry20160205/' ...
                                'GCAVisuals/Descriptor/Filopodia/Validation/' ...
                                'filoLength/ForMainMovie_Measurement_Movie'];
                            
                            file{2} =  [' ' projListC{iTraj,1} filesep imageFolder filesep '001.png '];
                            
                            
                            file{3} =  [' ' projListC{iTraj,1} filesep imageFolder filesep '061.png '];
                            
                            file{4} = [' ' projListC{iTraj,1} filesep imageFolder filesep '118.png '];
                            
                           
                        
                       % Make Overlay All 
                        arrayfun(@(i) scatter(y(currentTraj(i),1),y(currentTraj(i),2),50,cmapTime(i,:),'filled'),1:length(currentTraj));
                        arrayfun(@(i) plot(y(currentTraj(i):currentTraj(i+1),1),y(currentTraj(i):currentTraj(i+1),2),... 
                            'color',cmapTime(i,:)),1:length(currentTraj)-1);
                        % plot the start and end in bold 
%                         scatter(y(currentTraj(1),1),y(currentTraj(1),2),100,cmapTime(1,:),filled,'MarkerEdgeColor','w'); 
%                         
%                         scatter(y(currentTraj(end),1),y(currentTraj(end),2),100,cmapTime(end,:),filled,'MarkerEdgeColor','w'); 
                          %end
                        xlabel('MDS1');
                        ylabel('MDS2');
                        axis(axisLims);
                        
                        file{1} = [' ' perTrajDir filesep  'All.png ']; 

                        saveas(gcf,[perTrajDir filesep  'All.png']);
                        saveas(gcf,[perTrajDir filesep 'All.fig']); 
                        
                        
                        
                        filenamesStr = horzcat(file{:}); 
                        
                        montageDir  = [trajDir filesep 'Montages' ]; 
                        if ~isdir(montageDir) 
                          mkdir(montageDir); 
                        end 
                        nameSave = strrep(name,' ' , '_'); 
                        outC = [montageDir filesep 'Traj' nameSave '.png'];
                        
                        cmd = [' montage -tile ' 'x1' ' -geometry +5+5+0+0 -background "rgb(255,255,255)"' filenamesStr ...
                        '  -compress lzw ' outC];
                        system(cmd); 
                        
                    else % per group color 
                        scatter(y(currentTraj,1),y(currentTraj,2),50,[0.5,0.5,0.5],'filled','MarkerEdgeColor','k');
                        plot(y(currentTraj,1),y(currentTraj,2),'color',cmapCurrent(iTraj,:));
                        
                          %end
                        xlabel('MDS1');
                        ylabel('MDS2');
                        axis(axisLims);

                        
                    end
                    
                   
               
                
                close gcf
                else
                    display(['The Entire Movie' projListC{iTraj} 'was removed: Check Measurements']); 
                end 
                end 
               
            end
            %scatter(y(groupingTraj==i,1),y(:,2),50,'k','filled');
            
            % quick fix is htat all the trajectories should be 119
            % frames
            % overlay the trajectories colored by group (use the
            % shade colors shades)
            %% Optional- plot individual trajectories per group
            % Q and D for now : just need grouping trajectory (a list
            % of values marking the cell ID to which each frame belongs
            % and for iGroup = 1:nGroups
            %                   groupingCell = toPlot.info.perFrameGroupID; % rx1 vector of group IDs for each frame.
            %                   groupingTraj = toPlot.info.perFrameCellID; % rx1 vector of cell IDs for each frame in the compiled group set.
            %                   % where r is the total number of frames (observations)
            %                   % fed into the MDS
            %                   nTrajs = length(unique(groupingTraj(groupingCell==iGroup)));
            %                   for iTraj = 1:nTrajs
            %                       scatter(y(groupingTraj == iTraj  & groupingCell == iGroup,:), y(groupingTraj == iTraj & groupingCell == iGroup));
            %                   end
            
            
            %
            
            
%             
%             saveas(gcf,[ip.Results.OutputDirectory filesep 'MDS_PerFrame_2DScatter.fig']);
%             saveas(gcf,[ip.Results.OutputDirectory filesep 'MDS_PerFrame_2DScatter.png']);
%             saveas(gcf,[ip.Results.OutputDirectory filesep 'MDS_PerFrame_2DScatter.eps'],'psc2');
            
        else % plot by group
            
            arrayfun(@(x) scatter(y(grouping==x,1),...
                y(grouping==x,2),100,toPlot.info.color{x},'filled','MarkerEdgeColor',toPlot.info.color{x}),1:nGroups);
            hold on
            if ip.Results.plotCentroids
             
                %                       valuesCent = y(nCells + iGroup,:);
                arrayfun(@(x) scatter(y(nCells+x,1),y(nCells+x,2),200,toPlot.info.color{x},'+'),1:nGroups);
                hold on 
                             
            end 
            
            xlabel('MDS1');
            ylabel('MDS2');
            axis(axisLims);
            
            saveas(gcf,[cOutMDS filesep 'MDS_2DScatterByGroupColor.fig']);
            saveas(gcf,[cOutMDS filesep 'MDS_2DScatterByGroupColor.eps'],'psc2');
            saveas(gcf,[cOutMDS filesep 'MDS_2DScatterByGroupColor.png']);
        end % ip.Results.perFrame
        close gcf
           
        %% 
        if ip.Results.plotCentroids
            setAxis('on')
                hold on 
                
                for iNeurite = 1:size(IDsClosest,1)
                    img = imread([IDsClosest{iNeurite,1} filesep 'GrowthConeAnalyzer' filesep 'VisualizationOverlays' filesep 'Raw'...
                        filesep '119.png']);
                    valuesC = y(grouping == iNeurite,:);
                    image([valuesC(idxClosest{iNeurite},1)-0.4, valuesC(idxClosest{iNeurite},1)+0.4],[valuesC(idxClosest{iNeurite},2)-0.4,valuesC(idxClosest{iNeurite},2)+0.4],img);
                end
                
                
                %                       valuesCent = y(nCells + iGroup,:);
                arrayfun(@(x) scatter(y(nCells+x,1),y(nCells+x,2),200,toPlot.info.color{x},'+'),1:nGroups);
                
                
                for iNeurite = 1:size(IDsClosest,1)
                    img = imread([IDsFurthest{iNeurite,1} filesep 'GrowthConeAnalyzer' filesep 'VisualizationOverlays' filesep 'Raw'...
                        filesep '119.png']);
                    valuesC = y(grouping == iNeurite,:);
                    image([valuesC(idxFurthest{iNeurite},1)-0.4, valuesC(idxFurthest{iNeurite},1)+0.4],[valuesC(idxFurthest{iNeurite},2)-0.4,valuesC(idxFurthest{iNeurite},2)+0.4],img);
                end
                
                xlabel('MDS1');
                ylabel('MDS2');
                axis(axisLims);
            
            saveas(gcf,[cOutMDS filesep 'MDS_ViewLandscape.fig']);
            saveas(gcf,[cOutMDS filesep 'MDS_ViewLandscape.eps'],'psc2');
            saveas(gcf,[cOutMDS filesep 'MDS_ViewLandscape.png']);
              close gcf  
        end % plot Centroids
        
            
           
        %% plot control by cluster
        if ip.Results.cluster 
        cOutMDSClust = [cOutMDS filesep 'cluster'];
        if ~isdir(cOutMDSClust)
            mkdir(cOutMDSClust)
        end
        
        controlValues = y(grouping==1,:);
        %      groupingCluster = toPlot.info.groupingIdxClust;
        groupingControl = groupingCluster(grouping==1);
       
        nClusts = length(unique(groupingCluster(groupingCluster~=0)));
        setAxis('on')
        hold on
        arrayfun(@(x) scatter(controlValues(groupingControl==x,1),...
            controlValues(groupingControl==x,2),50, cMap(x,:),'filled'),1:nClusts);
        obsNames =  cellfun(@(x) strrep(x,'_',' '),obsNames,'uniformoutput',0);
        %      obsNames(test) = [];
        arrayfun(@(x) text(controlValues(x,1),controlValues(x,2),obsNames{x}),1:size(controlValues(:,1)));
        
        
        xlabel('MDS1');
        ylabel('MDS2');
        
        axis(axisLims);
        
        saveas(gcf,[cOutMDSClust filesep 'MDS_2DScatterByGroupClusterControl.fig']);
        saveas(gcf,[cOutMDSClust filesep 'MDS_2DScatterByGroupClusterControl.eps'],'psc2');
        saveas(gcf,[cOutMDSClust filesep 'MDS_2DScatterByGroupClusterControl.png']);
        close gcf
        
        %% Show Pics of Structures: Cluster 1 and Cluster 2
        name{1} = 'Low';
        name{2} = 'Middle';
        name{3} = 'High'; 
        for iClust = 1:nClusts
            setAxis('on')
            
            hold on
            
            %        scatter(controlValues(groupingControl==iClust,1),...
            %            controlValues(groupingControl==iClust,2),50, colorClust{iClust},'filled');
            
            arrayfun(@(x) scatter(controlValues(groupingControl==x,1),...
                controlValues(groupingControl==x,2),50,cMap(x,:),'filled'),1:nClusts);
            obsNames =  cellfun(@(x) strrep(x,'_',' '),obsNames,'uniformoutput',0);
            
            projListC = projList(groupingControl==iClust,1);
            valuesC = controlValues(groupingControl==iClust,:);
            for iNeurite = 1:size(projListC,1)
                img = imread([projListC{iNeurite,1} filesep 'GrowthConeAnalyzer' filesep 'VisualizationOverlays' filesep 'Raw'...
                    filesep '001.png']);
                image([valuesC(iNeurite,1)-0.4, valuesC(iNeurite,1)+0.4],[valuesC(iNeurite,2)-0.4,valuesC(iNeurite,2)+0.4],img);
                clear img
            end
            title([name{iClust} ' Outgrowth Control Structures Visualized']);
            
            xlabel('MDS1');
            ylabel('MDS2');
            
            axis(axisLims);
            saveas(gcf,[cOutMDSClust filesep 'MDS_2DVisualsCluster_' name{iClust} '.fig']);
            saveas(gcf,[cOutMDSClust filesep 'MDS_2DVisualsCluster_' name{iClust} '.eps'],'psc2');
            saveas(gcf,[cOutMDSClust filesep 'MDS_2DVisualsCluster_' name{iClust} '.png']);
            close gcf
            
        end
        
        %% plot all by cluster
        setAxis('on')
        hold on
        arrayfun(@(x) scatter(y(groupingCluster==x,1),...
            y(groupingCluster==x,2),100,cMap(x,:),'filled'),1:nClusts);
        
        
        %arrayfun(@(x) text(y(x,1),y(x,2),obsNames{x}),1:size(y(:,1)))
        xlabel('MDS1');
        ylabel('MDS2');
        
        axis(axisLims);
        
        saveas(gcf,[cOutMDSClust filesep 'MDS_2DScatterByGroupClusterAll.fig']);
        saveas(gcf,[cOutMDSClust filesep 'MDS_2DScatterByGroupClusterAll.eps'],'psc2');
        saveas(gcf,[cOutMDSClust filesep 'MDS_2DScatterByGroupClusterAll.png']);
        close gcf
        end 
        
         
        %% plot results by outgrowth scaleMap
        if ip.Results.plotByGroupElongation
            %
            MDSGroupDir = [cOutMDS filesep 'MDSPerGroup'];
            if ~isdir(MDSGroupDir);
                mkdir(MDSGroupDir);
            end
            
            setAxis('on');
           
              for iGroup = 1:nGroups
                
                
              
%                 if iGroup > 1
%                     scatter(controlValues(:,1),controlValues(:,2),50,'k','filled');
%                 end
%                 
                % get group values
                yC = y(grouping ==iGroup,:);
                
                idxCMapC = idxCMap(grouping == iGroup,:);
                
                % Plot Current group in color
                for iColor = 1:length(cmap)
                    if ~isempty(idxCMapC(:,1)==iColor)
                        scatter((yC(idxCMapC(:,1) == iColor,1)),...
                            (yC(idxCMapC(:,1) == iColor,2)),100,cmap(iColor,:),'filled','MarkerEdgeColor',[0 0 0]);
                        hold on
                    end
                end
                
                
                axis(axisLims);
                %axis([-8,8,-8,8]);
%                 title(['Group ' toPlot.info.names{iGroup}]);
                xlabel('MDS1');
                ylabel('MDS2');
                saveas(gcf,[MDSGroupDir filesep 'MDS_2DScatterColorByOutgrowthAll.fig']);
                saveas(gcf,[MDSGroupDir filesep 'MDS_2DScatterColorByOutgrowthAll.eps'],'psc2');
                saveas(gcf,[MDSGroupDir filesep 'MDS_2DScatterColorByOutgrowthAll.png']);        
              
              end 
            
            
            
            for iGroup = 1:nGroups
                setAxis('on');
                hold on
                
                
                
                
                
                
                if iGroup > 1
                    scatter(controlValues(:,1),controlValues(:,2),50,'k','filled');
                end
                
                % get group values
                yC = y(grouping ==iGroup,:);
                
                idxCMapC = idxCMap(grouping == iGroup,:);
                
                % Plot Current group in color
                for iColor = 1:length(cmap)
                    if ~isempty(idxCMapC(:,1)==iColor)
                        scatter((yC(idxCMapC(:,1) == iColor,1)),...
                            (yC(idxCMapC(:,1) == iColor,2)),100,cmap(iColor,:),'filled','MarkerEdgeColor',[0 0 0]);
                    end
                end
                % Plot the rest in gray
                
                
                %axis([min(y(:,1)),max(y(:,1)),min(y(:,2)),max(y(:,2))]);
                axis(axisLims);
                %axis([-8,8,-8,8]);
                title(['Group ' toPlot.info.names{iGroup}]);
                xlabel('MDS1');
                ylabel('MDS2');
                saveas(gcf,[MDSGroupDir filesep 'MDS_2DScatterColorByOutgrowthAll_' toPlot.info.names{iGroup} '.fig']);
                saveas(gcf,[MDSGroupDir filesep 'MDS_2DScatterColorByOutgrowthAll_' toPlot.info.names{iGroup} '.eps'],'psc2');
                saveas(gcf,[MDSGroupDir filesep 'MDS_2DScatterColorByOutgrowthAll_' toPlot.info.names{iGroup} '.png']);
                close gcf
            end % iGroup 
        end % plot by group elongation 
                %% plot results by group with picture 
                
        if ip.Results.plotByGroupPic
                MDSGroupDir = [cOutMDS filesep 'MDSPerGroupPic'];
                if ~isdir(MDSGroupDir); 
                    mkdir(MDSGroupDir)
                end 
                
                for iGroup = 2:nGroups
                      
                setAxis('on');
                hold on
                
           
                % Plot the Control by outgrowth first 
                yC = y(grouping ==1,:);
                
                idxCMapC = idxCMap(grouping == 1,:);
                
                % Plot Current group in color
                for iColor = 1:length(cmap)
                    if ~isempty(idxCMapC(:,1)==iColor)
                        scatter((yC(idxCMapC(:,1) == iColor,1)),...
                            (yC(idxCMapC(:,1) == iColor,2)),100,cmap(iColor,:),'filled','MarkerEdgeColor',[0 0 0]);
                    end
                end
                hold on 
                
                
%                     yC = y(grouping ==iGroup,:);
%                     idxCMapC = idxCMap(grouping == iGroup,:);
%                     for iColor = 1:length(cmap)
%                         if ~isempty(idxCMapC(:,1)==iColor)
%                             % get the coords
%                             x = yC(idxCMap(:,1)==iColor); 
%                             y = yC(idxCMap(:,1)==iColor); 
%                             img = imread([projListC{iNeurite,1} filesep 'GrowthConeAnalyzer' filesep 'VisualizationOverlays' filesep 'Raw'...
%                             filesep '001.png']);
%                             image([x-0.4, x+0.4],[y-0.4,y+0.4],img);
%                             scatter(
% %                             scatter((yC(idxCMapC(:,1) == iColor,1)),...
% %                                 (yC(idxCMapC(:,1) == iColor,2)),100,cmap(iColor,:),'filled','MarkerEdgeColor',[0 0 0]);
%                         end
%                     end
                    projListC = toPlot.info.projList{iGroup}; 
                  
                    valuesC = y(grouping==iGroup,:); 
                    for iNeurite = 1:size(projListC,1)
                        img = imread([projListC{iNeurite,1} filesep 'GrowthConeAnalyzer' filesep 'VisualizationOverlays' filesep 'Raw'...
                            filesep '001.png']);
                        image([valuesC(iNeurite,1)-0.4, valuesC(iNeurite,1)+0.4],[valuesC(iNeurite,2)-0.4,valuesC(iNeurite,2)+0.4],img);
                        clear img
                    end
                    %axis([min(y(:,1)),max(y(:,1)),min(y(:,2)),max(y(:,2))]);
                    axis(axisLims);
                    %axis([-8,8,-8,8]);
                    title(['Group ' toPlot.info.names{iGroup}]);
                    xlabel('MDS1');
                    ylabel('MDS2');
                    
                    saveas(gcf,[MDSGroupDir filesep 'MDS_2DScatterImagePerGroup_' toPlot.info.names{iGroup} '.fig']);
                    saveas(gcf,[MDSGroupDir filesep 'MDS_2DScatterImagePerGroup_' toPlot.info.names{iGroup} '.eps'],'psc2');
                    saveas(gcf,[MDSGroupDir filesep 'MDS_2DScatterImagePerGroup_' toPlot.info.names{iGroup} '.png']);
                    
                    close gcf
                    
                end % iGroup 
        end % plotByGroupPic    

%         
        %% plot by group color 
          if ip.Results.plotByGroupColorID    
              MDSGroupDir = [cOutMDS filesep 'MDSPerGroupColorByGroup'];
              if ~isdir(MDSGroupDir);
                  mkdir(MDSGroupDir)
              end
              
              for iGroup = 2:nGroups
                  
                  setAxis('on',0.75,24);  
                  
                  % Plot the Control by in black first 
                  yC = y(grouping ==1,:);
                  
                  
                  if ip.Results.plotCentroids
                      valuesCent = y(nCells +1,:);
                      scatter(valuesCent(:,1),valuesCent(:,2),200,'k','+'); 
                      %% SECTION TITLE
                      % DESCRIPTIVE TEXT
                      
                      
                  end 
                  
                  scatter(yC(:,1),yC(:,2),100,'k','filled'); 
                  hold on 
                  
                  valuesC = y(grouping==iGroup,:);
                  scatter(valuesC(:,1),valuesC(:,2),300,toPlot.info.color{iGroup},'filled');
                  
                  
                  if ip.Results.plotCentroids
                      valuesCent = y(nCells + iGroup,:);
                      scatter(valuesCent(:,1),valuesCent(:,2),200,toPlot.info.color{iGroup},'+'); 
                      
                     
                  end 
                
                  axis(axisLims);
                  %axis([-8,8,-8,8]);
                  title(['Group ' toPlot.info.names{iGroup}]);
                  xlabel('MDS1');
                  ylabel('MDS2');
                 
                  saveas(gcf,[MDSGroupDir filesep 'MDS_2DScatterImagePerGroup_' toPlot.info.names{iGroup} '.fig']);
                  saveas(gcf,[MDSGroupDir filesep 'MDS_2DScatterImagePerGroup_' toPlot.info.names{iGroup} '.eps'],'psc2');
                  saveas(gcf,[MDSGroupDir filesep 'MDS_2DScatterImagePerGroup_' toPlot.info.names{iGroup} '.png']);
                  
                  close gcf
                  
              end % iGroup
              
              
              
              if ip.Results.cluster
                  clustDirElongNorm = [ cOutMDS filesep 'cluster' filesep 'CompareKDToElongNorm'];
                  if ~isdir(clustDirElongNorm)
                      mkdir(clustDirElongNorm);
                  end
                  
                  for iGroup = 2:nGroups
                      setAxis('on');
                      
                      % Plot the Control by in black first
                      yC = y(grouping ==1 & groupingCluster ==2 ,:);
                      
                      scatter(yC(:,1),yC(:,2),50,'k','filled');
                      hold on
                      
                      valuesC = y(grouping==iGroup,:);
                      scatter(valuesC(:,1),valuesC(:,2),100,toPlot.info.color{iGroup},'filled');
                      
                      axis(axisLims);
                      %axis([-8,8,-8,8]);
                      title(['Group ' toPlot.info.names{iGroup}]);
                      xlabel('MDS1');
                      ylabel('MDS2');
                      
                      saveas(gcf,[ clustDirElongNorm filesep 'CompareNormElongControl_vs_' toPlot.info.names{iGroup} '.fig']);
                      saveas(gcf,[ clustDirElongNorm filesep 'CompareNormElongControl_vs_' toPlot.info.names{iGroup} '.eps'],'psc2');
                      saveas(gcf,[ clustDirElongNorm filesep 'CompareNormElongControl_vs_' toPlot.info.names{iGroup} '.png']);
                      
                      close all
                      
                      %                   saveas(gcf,[MDSGroupDir filesep 'MDS_2DScatterImagePerGroup_' toPlot.info.names{iGroup} '.fig']);
                      %                   saveas(gcf,[MDSGroupDir filesep 'MDS_2DScatterImagePerGroup_' toPlot.info.names{iGroup} '.eps'],'psc2');
                      %                   saveas(gcf,[MDSGroupDir filesep 'MDS_2DScatterImagePerGroup_' toPlot.info.names{iGroup} '.png']);
                      
                  end
                  
                  clustDirElongLow = [ cOutMDS filesep 'cluster' filesep 'CompareKDToElongLow'];
                  if ~isdir(clustDirElongLow)
                      mkdir(clustDirElongLow);
                  end
                  
                  for iGroup = 2:nGroups
                      setAxis('on');
                      
                      % Plot the Control by in black first
                      yC = y(grouping ==1 & groupingCluster ==1,:);
                      
                      scatter(yC(:,1),yC(:,2),50,'k','filled');
                      hold on
                      
                      valuesC = y(grouping==iGroup,:);
                      scatter(valuesC(:,1),valuesC(:,2),100,toPlot.info.color{iGroup},'filled');
                      
                      axis(axisLims);
                      %axis([-8,8,-8,8]);
                      title(['Group ' toPlot.info.names{iGroup}]);
                      xlabel('MDS1');
                      ylabel('MDS2');
                      
                      saveas(gcf,[ clustDirElongLow filesep 'CompareLowElongControl_vs_' toPlot.info.names{iGroup} '.fig']);
                      saveas(gcf,[ clustDirElongLow filesep 'CompareLowElongControl_vs_' toPlot.info.names{iGroup} '.eps'],'psc2');
                      saveas(gcf,[ clustDirElongLow filesep 'CompareLowElongControl_vs_' toPlot.info.names{iGroup} '.png']);
                      
                      close all
                      
                      %                   saveas(gcf,[MDSGroupDir filesep 'MDS_2DScatterImagePerGroup_' toPlot.info.names{iGroup} '.fig']);
                      %                   saveas(gcf,[MDSGroupDir filesep 'MDS_2DScatterImagePerGroup_' toPlot.info.names{iGroup} '.eps'],'psc2');
                      %                   saveas(gcf,[MDSGroupDir filesep 'MDS_2DScatterImagePerGroup_' toPlot.info.names{iGroup} '.png']);
                      
                  end
                  
                  
         
              end
              
              
              
              
              
              
              
              
              
              
              
              
              
              
          end % ip.Results.plotByGroupColorID   
              
    end      %% for iCriterion 
              %% plot for each group
      
end % ip.Results.MDS


%% Perform the PCA
if ip.Results.PCA

    PCADir = [ip.Results.OutputDirectory filesep 'PCA']; 
    if ~isdir(PCADir)
        mkdir(PCADir)
    end 

    [coef,scores,latent,ts,exp] = pca(dataFinal);
    coef = double(coef);
    scores = double(scores);
    save([PCADir filesep 'PCAData.mat'],'coef','scores','latent','ts','exp');
    save([PCADir filesep 'PCAZScoresFinal.mat'],'dataFinal');
    
    
    

    %% Plot variance explained.
    setAxis
    
    cumExp = cumsum(latent)./sum(latent);
    PC = 1:length(cumExp);
    scatter(PC,cumExp,20,'filled','k');
    ylabel('Variance Explained');
    xlabel('PC Number');
    forLine = PC(find(cumExp>.95,1,'first'));
    line([forLine,forLine],[0:1],'color','k');
    
    saveas(gcf,[PCADir filesep 'PercentVarianceExplained.fig']);
    saveas(gcf,[PCADir filesep 'PercentVarianceExplained.eps'],'psc2');
    saveas(gcf,[PCADir filesep 'PercentVarianceExplained.png']);
    
    %% 2D Plots
    for iPC = 1:6;
        
        cDir = ([PCADir filesep 'PC' num2str(iPC) 'vs' 'PC' num2str(iPC+1) ] ) ;
        if ~isdir(cDir)
            mkdir(cDir);
        end
        %% Biplot
        
        fsFigure(0.75)
        biplot(coef(:,iPC:iPC+1),'varlabels',varNamesString,'obslabels',obsNames);
        xlabel({['PC' num2str(iPC)] ; ['Percent Variance Explained ' num2str(exp(iPC),3) '%']});
        ylabel({['PC' num2str(iPC+1)] ; ['Percent Variance Explained ' num2str(exp(iPC+1),3) '%']});
        
        saveas(gcf,[cDir filesep 'PCABiPlotPC' num2str(iPC) 'vs' num2str(iPC+1) '.fig']);
        saveas(gcf,[cDir filesep 'PCABiPlotPC' num2str(iPC) 'vs' num2str(iPC+1) '.eps'],'psc2');
        saveas(gcf,[cDir filesep 'PCABiPlotPC' num2str(iPC) 'vs' num2str(iPC+1) '.png']);
        
        %% Scores (Per Neurite) Overlay ColorCoded By Experimental Condition
        scaleFact1 = max(max(abs(scores(:,iPC:iPC+1)))) ;
        coefs = coef(:,iPC:iPC+1);
        maxCoefLen = sqrt(max(sum(coefs.^2,2)));
        
        hold on
        arrayfun(@(x) scatter((scores(grouping==x,iPC)./scaleFact1).*maxCoefLen,...
            (scores(grouping==x,iPC+1)./scaleFact1).*maxCoefLen,50,toPlot.info.color{x},'filled'),1:nGroups);
        
        grid('off')
        saveas(gcf,[cDir filesep 'PCABiPlotPC' num2str(iPC) 'vs' num2str(iPC+1) '_ScoresByGroup.fig']);
        saveas(gcf,[cDir filesep 'PCABiPlotPC' num2str(iPC) 'vs' num2str(iPC+1) '_ScoresByGroup.eps'],'psc2');
        saveas(gcf,[cDir filesep 'PCABiPlotPC' num2str(iPC) 'vs' num2str(iPC+1) '_ScoresByGroup.png']);
        
        %% 
         arrayfun(@(x) scatter((scores(grouping==x,iPC)./scaleFact1).*maxCoefLen,...
            (scores(grouping==x,iPC+1)./scaleFact1).*maxCoefLen,50,toPlot.info.color{x},'filled'),1:nGroups);
         close gcf
        %%
        %     % Plot Events Selected
        %     for iColor = 1:length(cmap)
        %         if ~isempty(idxCMap(:,1)==iColor)
        %             scatter((scores(idxCMap(:,1) == iColor,1)./scaleFact1).*maxCoefLen,...
        %                 (scores(idxCMap(:,1) == iColor,2)./scaleFact1).*maxCoefLen,50,'MarkerEdgeColor',cmap(iColor,:));
        %         end
        %     end
        %
        PCAImByGroupDir = [cDir filesep 'StructByGroup'] ;
        if ~isdir(PCAImByGroupDir)
            mkdir(PCAImByGroupDir);
        end
        for iGroup = 1:nGroups
            
            
            fsFigure(0.75)
            biplot(coef(:,iPC:iPC+1),'varlabels',varNamesString,'obslabels',obsNames);
            xlabel({['PC' num2str(iPC)] ; ['Percent Variance Explained ' num2str(exp(iPC),3) '%']});
            ylabel({['PC' num2str(iPC+1)] ; ['Percent Variance Explained ' num2str(exp(iPC+1),3) '%']});
            grid ('off')
            hold on
            scoresC = scores(grouping ==iGroup,:)./scaleFact1.*maxCoefLen;
            projListC = toPlot.info.projList{iGroup};
            for iNeurite = 1:size(scoresC,1);
                xC = scoresC(iNeurite,iPC);
                yC = scoresC(iNeurite,iPC+1);
                
                img = imread([projListC{iNeurite,1} filesep 'GrowthConeAnalyzer' filesep 'VisualizationOverlays' filesep 'Raw'...
                    filesep '001.png']);
                image([xC-0.04, xC+0.04],[yC-0.04,yC+0.04],img)
            end % iNeurite
            
            
              title(['Group ' toPlot.info.names{iGroup}]);
                  
                    
                    saveas(gcf,[  PCAImByGroupDir filesep 'PCABiPlotPC' num2str(iPC) 'vs' num2str(iPC+1) '_' toPlot.info.names{iGroup} '.fig']);
                    saveas(gcf,[  PCAImByGroupDir filesep 'PCABiPlotPC' num2str(iPC) 'vs' num2str(iPC+1) '_' toPlot.info.names{iGroup} '.eps'],'psc2');
                    saveas(gcf,[  PCAImByGroupDir filesep 'PCABiPlotPC' num2str(iPC) 'vs' num2str(iPC+1) '_' toPlot.info.names{iGroup} '.png']);
   close gcf
        end % iGroup
         
        
       
        %% Scores (Per Neurite) Overlay ColorCoded By Outgrowth
        
        fsFigure(0.75)
        biplot(coef(:,iPC:iPC+1),'varlabels',varNamesString,'obslabels',obsNames);
        
        hold on
        % Plot Events Selected
        for iColor = 1:length(cmap)
            if ~isempty(idxCMap(:,1)==iColor)
                scatter((scores(idxCMap(:,1) == iColor,iPC)./scaleFact1).*maxCoefLen,...
                    (scores(idxCMap(:,1) == iColor,iPC+1)./scaleFact1).*maxCoefLen,50,cmap(iColor,:),'filled','MarkerEdgeColor',[0 0 0]);
            end
        end
        
        xlabel({['PC' num2str(iPC)] ; ['Percent Variance Explained ' num2str(exp(iPC),3) '%']});
        ylabel({['PC' num2str(iPC+1)] ; ['Percent Variance Explained ' num2str(exp(iPC+1),3) '%']});
        
        saveas(gcf,[cDir filesep 'PCABiPlotPC' num2str(iPC) 'vs' num2str(iPC+1) '_ScoresByNetOutgrowth10Min.fig']);
        saveas(gcf,[cDir filesep 'PCABiPlotPC' num2str(iPC) 'vs' num2str(iPC+1) '_ScoresByNetOutgrowth10Min.eps'],'psc2');
        saveas(gcf,[cDir filesep 'PCABiPlotPC' num2str(iPC) 'vs' num2str(iPC+1) '_ScoresByNetOutgrowth10Min.png']);
        
        hold on
        obsNames =  cellfun(@(x) strrep(x,'_',' '),obsNames,'uniformoutput',0);
        
        arrayfun(@(x) text((scores(x,iPC)./scaleFact1).*maxCoefLen,(scores(x,iPC+1)./scaleFact1).*maxCoefLen,obsNames{x}),1:length(scores(:,1)));
        
        saveas(gcf,[cDir filesep 'PCABiPlotPC' num2str(iPC) 'vs' num2str(iPC+1) '_ScoresByNetOutgrowth10MinWithText.fig']);
        saveas(gcf,[cDir filesep 'PCABiPlotPC' num2str(iPC) 'vs' num2str(iPC+1) '_ScoresByNetOutgrowth10MinWithText.eps'],'psc2');
        saveas(gcf,[cDir filesep 'PCABiPlotPC' num2str(iPC) 'vs' num2str(iPC+1) '_ScoresByNetOutgrowth10MinWithText.png']);
        close gcf
    end % for iPC
    
    %% 3D data plots
    %% Biplot
    fsFigure(0.75)
    biplot(coef(:,1:3),'varlabels',varNamesString);
    xlabel({'PC1' ; [ 'Percent Variance Explained ' num2str(exp(1),3) '%']});
    ylabel({'PC2' ; [ 'Percent Variance Explained ' num2str(exp(2),3) '%']});
    zlabel({'PC3' ; [ 'Percent Variance Explained ' num2str(exp(3),3) '%']});
    axis([-0.5 0.5 -0.5 0.5 -0.5 0.5]);
    view(45,45);
    
    saveas(gcf,[PCADir filesep 'PCABiPlotFirst3PCs.fig']);
    saveas(gcf,[PCADir 'PCABiPlotFirst3PCs.eps'],'psc2');
    saveas(gcf,[PCADir filesep 'PCABiPlotFirst3PCs.png']);
    
    fsFigure(0.75)
    %% Color by Group Condition
    
    biplot(coef(:,1:3));
    scaleFact1 = max(max(abs(scores(:,1:3)))) ;
    coefs = coef(:,1:3);
    maxCoefLen = sqrt(max(sum(coefs.^2,2)));
    
    hold on
    arrayfun(@(x) scatter3((scores(grouping==x,1)./scaleFact1).*maxCoefLen,...
        (scores(grouping==x,2)./scaleFact1).*maxCoefLen,(scores(grouping==x,3)./scaleFact1).*maxCoefLen,...
        50,toPlot.info.color{x},'filled'),1:nGroups);
    
    xlabel({'PC1' ; ['Percent Variance Explained ' num2str(exp(1),3) '%']});
    ylabel({'PC2' ; ['Percent Variance Explained ' num2str(exp(2),3) '%']});
    zlabel({'PC3' ; ['Percent Variance Explained ' num2str(exp(3),3) '%']});
    
    axis([-0.5 0.5 -0.5 0.5 -0.5 0.5]);
    view(45,45);
    saveas(gcf,[PCADir  filesep 'PCABiPlotFirst3PCs_scatter.fig']);
    saveas(gcf,[PCADir filesep 'PCABiPlotFirst3PCs_scatter.eps'],'psc2');
    saveas(gcf,[PCADir filesep 'PCABiPlotFirst3PCs_scatter.png']);
    close gcf
    
    %%
    biplot(coef(:,1:3));
    hold on
    for iColor = 1:length(cmap)
        if ~isempty(idxCMap(:,1)==iColor)
            scatter3((scores(idxCMap(:,1) == iColor,1)./scaleFact1).*maxCoefLen,...
                (scores(idxCMap(:,1) == iColor,2)./scaleFact1).*maxCoefLen, ...
                (scores(idxCMap(:,1) == iColor,3)./scaleFact1).*maxCoefLen,...
                50,cmap(iColor,:),'filled','MarkerEdgeColor',[0 0 0]);
        end
    end
    xlabel({'PC1' ; ['Percent Variance Explained ' num2str(exp(1),3) '%']});
    ylabel({'PC2' ; ['Percent Variance Explained ' num2str(exp(2),3) '%']});
    zlabel({'PC3' ; ['Percent Variance Explained ' num2str(exp(3),3) '%']});
    
    view(45,45);
    axis([-0.5 0.5 -0.5 0.5 -0.5 0.5])
    saveas(gcf,[PCADir  filesep 'PCABiPlotFirst3PCs_scatterColorByGrowth.fig']);
    saveas(gcf,[PCADir filesep 'PCABiPlotFirst3PCs_scatterColorByGrowth.eps'],'psc2');
    saveas(gcf,[PCADir filesep 'PCABiPlotFirst3PCs_scatterColorByGrowth.png']);
    
    %% Plot each group by color
    if ip.Results.plotByGroup
        if ~isempty(ip.Results.OutputDirectory);
            outDirByGroup = [ip.Results.OutputDirectory filesep 'perGroupScatters'];
            if ~isdir(outDirByGroup)
                mkdir(outDirByGroup)
            end
            
        end
        
        for iGroup = 1:nGroups
            % subplot(2,4,iGroup);  % too small
            fsFigure(0.75)
            biplot(coef(:,1:3));
            if strcmpi(toPlot.info.names{iGroup}, 'KDNo');
                forTitle = 'Control';
            else
                forTitle = toPlot.info.names{iGroup};
            end
            title(forTitle,'FontName','Arial','FontSize',22);
            hold on
            
            
            for iColor = 1:length(cmap)
                if ~isempty(idxCMap(:,1)==iColor)
                    
                    scatter3((scores(idxCMap(:,1) == iColor & grouping ==iGroup,1)./scaleFact1).*maxCoefLen,...
                        (scores(idxCMap(:,1) == iColor & grouping ==iGroup,2)./scaleFact1).*maxCoefLen, ...
                        (scores(idxCMap(:,1) == iColor & grouping == iGroup,3)./scaleFact1).*maxCoefLen,...
                        50,cmap(iColor,:),'filled','MarkerEdgeColor',[0 0 0]);
                end
            end % iColor
            view(45,45);
            axis([-0.5 0.5 -0.5 0.5 -0.5 0.5])
            xlabel({'PC1' ; ['Percent Variance Explained ' num2str(exp(1),3) '%']});
            ylabel({'PC2' ; ['Percent Variance Explained ' num2str(exp(2),3) '%']});
            zlabel({'PC3' ; ['Percent Variance Explained ' num2str(exp(3),3) '%']});
            
            saveas(gcf,[outDirByGroup filesep num2str(iGroup,'%02d') '_' toPlot.info.names{iGroup} '.fig']);
            saveas(gcf,[outDirByGroup filesep num2str(iGroup,'%02d') '_' toPlot.info.names{iGroup} '.png']);
            saveas(gcf,[outDirByGroup filesep num2str(iGroup,'%02d') '_' toPlot.info.names{iGroup} '.eps'],'psc2');
            close gcf
        end % iGroup
    end
    
end % if ip.Result.PCA
%% Heterogeneity Plots : average distance of cells from centroid

if ip.Results.DistMetrics
    
    distanceDir = [ip.Results.OutputDirectory filesep 'DistanceMetrics']; 
    if ~isdir(distanceDir) 
        mkdir(distanceDir); 
    end 
    
    if ip.Results.cluster
        setAxis('on'); 
        % total scatter 

        dataControl = dataFinal(grouping ==1,:); 
      
               
        % Cluster By elongation a response variable. 
        % Calculate the euclidean distance from each point in the
        % descriptor space and the centroid of this cluster. 
        
        % get each cluster
        numClust = length(unique(groupingCluster(groupingCluster~=0))); 
        data{1} = dataControl; 
        
        dataByClust = arrayfun(@(x) dataFinal(groupingCluster==x & grouping ==1,:),1:numClust,'uniformoutput',0);
        projListByClust = arrayfun(@(x) projList(groupingCluster==x & grouping ==1,:),1:numClust,'uniformoutput',0); 
        save([distanceDir filesep 'projListByClust.mat'],'projListByClust');
      
        
        data = [data,dataByClust];
        %[meanValue,scatterValues] = arrayfun(@(x) calcMeanDistToCent(dataMatAllMeas(groupingCluster == x),:)',1:nClust,'uniformoutput',0) ;
        
        % get the confidence intervals of each cluster mean scatter (output a
        % cell with CIs cluster 1:nClust
        CIs = cellfun(@(x) bootci(2000,@calcMeanDistToCent,x'),data,'uniformoutput',0);
        hold on 
        
        
        [meanValue, scatterValues] = cellfun(@(x) calcMeanDistToCent(x'),data,'uniformoutput',0);
        %save([distanceDir filesep 'scatterValues.mat'],'scatterValues'); 
        
        scatterByClust =  arrayfun(@(x) [projListByClust{x} num2cell(scatterValues{x+1})],1:numClust,'uniformoutput',0); 
        save([distanceDir filesep 'scatterByClust.mat'],'scatterByClust'); 
        title('Mean Of Scatter','FontSize',10); 
       
        test = reformatDataCell(scatterValues);
        h =  notBoxPlot(test);
        hold on
        cMap = [[0,0,0]; cMap];
        arrayfun(@(i) errorbar(i,meanValue{i},meanValue{i} - CIs{i}(1) , CIs{i}(2) - meanValue{i},'color',cMap(i,:)),1:numel(CIs)); 
        arrayfun(@(i) set(h(i).semPtch,'faceColor','w'),1:numel(scatterValues));
        arrayfun(@(i) set(h(i).sdPtch,'faceColor','w'),1:numel(scatterValues));
        arrayfun(@(i) set(h(i).data,'markerFaceColor',cMap(i,:)),1:numel(scatterValues));
        arrayfun(@(i) set(h(i).mu,'color',cMap(i,:)),1:numel(scatterValues));
        %             arrayfun(@(i) set(h(i).data,'markerEdgeColor','w'),1:numel(toPlot.info.names));
        %             arrayfun(@(i) set(h(i).mu,'color',toPlot.info.color{i}),1:numel(toPlot.info.names));
        ylabel('Euclidean Distance From Centroid');
        
        
        names{1} = 'All'; 
        names{2} = 'Low '  ;     
        names{3} = 'Mid';
        names{4} = 'High '; 
        
        set(gca,'XTick',1:numel(names));
        set(gca,'XTickLabel',names,'FontSize',10);
       
        
        saveas(gcf,[distanceDir filesep  'ScatterMetricByElongCluster.png'] ); 
        saveas(gcf, [distanceDir filesep 'ScatterMetricByElongCluster.eps'],'psc2'); 
        saveas(gcf, [distanceDir filesep 'ScatterMetricByElongCluster.fig']); 
        
        close gcf
        
        %% for now put median in a separate section : eventually consolidate 
    
        
          setAxis('on'); 
 
        
        % get the confidence intervals of each cluster median scatter (output a
        % cell with CIs cluster 1:nClust
        CIs = cellfun(@(x) bootci(2000,@calcMedianDistToCent,x'),data,'uniformoutput',0);
        hold on 
        
        
        [medianValue, scatterValues] = cellfun(@(x) calcMedianDistToCent(x'),data,'uniformoutput',0);
      
       
        test = reformatDataCell(scatterValues);
        h =  notBoxPlot(test);
        hold on
%         cMap = [[0,0,0]; cMap];
        arrayfun(@(x) scatter(x,medianValue{x},100,cMap(x,:),'x'),1:numel(CIs)); 
        arrayfun(@(i) errorbar(i,medianValue{i},medianValue{i} - CIs{i}(1) , CIs{i}(2) - medianValue{i},'color',cMap(i,:)),1:numel(CIs)); 
        arrayfun(@(i) set(h(i).semPtch,'faceColor','w'),1:numel(scatterValues));
        arrayfun(@(i) set(h(i).sdPtch,'faceColor','w'),1:numel(scatterValues));
        arrayfun(@(i) set(h(i).data,'markerFaceColor',cMap(i,:)),1:numel(scatterValues));
        arrayfun(@(i) set(h(i).mu,'color','w'),1:numel(scatterValues));
        %             arrayfun(@(i) set(h(i).data,'markerEdgeColor','w'),1:numel(toPlot.info.names));
        %             arrayfun(@(i) set(h(i).mu,'color',toPlot.info.color{i}),1:numel(toPlot.info.names));
        ylabel('Euclidean Distance From Centroid');
        names{1} = 'All'; 
        names{2} = 'Low '  ;     
        names{3} = 'Mid';
        names{4} = 'High '; 
        set(gca,'XTick',1:numel(names));
        set(gca,'XTickLabel',names,'FontSize',10);
        title('Median Of Scatter','FontSize',10); 
%         
        saveas(gcf,[distanceDir filesep  'ScatterMetricMedianByElongCluster.png'] ); 
        saveas(gcf, [distanceDir filesep 'ScatterMetricMedianByElongCluster.eps'],'psc2'); 
        saveas(gcf, [distanceDir filesep 'ScatterMetricMedianByElongCluster.fig']);
        
        close gcf
        
        
    end
       
    
    
    
    %% Discrimination Metrics
    nCond = numel(dataMatAllMeas);
    controlCompare{1} = dataFinalAll{1};
    if ip.Results.cluster
        groupingControl = groupingCluster(grouping==1);
        controlCompare{2} = dataFinalAll{1}(groupingControl==1,:);
        controlCompare{3} = dataFinalAll{1}(groupingControl==2,:);
        add{2} = 'LowElongRate';
        add{3} = 'HighElongRate';
    end
   
    
    
    add{1} = 'All'; 
   
    
    for iCompare = 1:numel(controlCompare)
        
        %[DB,Dunn,SC,c1,c2] = arrayfun(@(x) whDiscriminationMeasures(dataMatAllMeas{x}',dataMatAllMeas{1}'),1:nCond);
        [DB,pValues] = arrayfun(@(x) bootDiscrimMetrics(dataFinalAll{x}',controlCompare{iCompare}','nrep',5000),1:nCond,'uniformoutput',0);
        %[DB,pValues] = arrayfun(@(x) bootDiscrimMetrics(dataMatAllMeas{x}',dataMatAllMeas{1}','nrep',5000),1:nCond,'uniformoutput',0);
        DB = vertcat(DB{:})';
        pValues = vertcat(pValues{:});
        
        setAxis('on');
        [DBSort,idxSort] = sort(DB);
        DBSort = DBSort(2:end);
        pValues = pValues(idxSort,:);
        pValues = pValues(2:end,:);
        
        
        colorsGroup = toPlot.info.color;
        colorsGroup = colorsGroup(idxSort);
        colorsGroup = colorsGroup(2:end);
        
        %     arrayfun(@(x) scatter(x,DBSort(x),50,colorsGroup(x,:),'filled'),1:length(DBSort));
        %     errorbar(@(x) errorbar(x,DBSort(x),CIs(x,1),CIs(x,2),'color',colorsGroup(x,:)),1:length(DBSort));
        % arrayfun(@(x)
        
        
        setAxis('on');
        h = bar([DBSort(1:end); DBSort(1:end)]); % do two for now and crop one out
        hold on
        for x = 1:length(h)
            h(x).FaceColor = colorsGroup{x};
            %errorbar(x,DB(x),CIs(x,1),CIs(x,2));
        end
        
        axis([1.5,2.5,0,max(DBSort(2:end))+0.2])
        
        namesGroup = toPlot.info.names;
        namesSort = namesGroup(idxSort);
        
        ylabel({'Separation Statistic-'; 'Inverse Davies Bouldin Index'});
        set(gca,'XTick',1:numel(namesSort)-1);
        set(gca,'XTickLabel',' ');
        
        saveas(gcf , [distanceDir filesep  'InverseDBIndexWithGroup_' add{iCompare} '.png']);
        saveas(gcf , [distanceDir filesep  'InverseDBIndexWithGroup_' add{iCompare} '.eps'],'psc2');
        saveas(gcf , [distanceDir filesep  'InverseDBIndexWithGroup_' add{iCompare} '.fig']);
        save([distanceDir filesep 'pValues_' add{iCompare}],'pValues');
        close gcf
    end 

   
    %% 
    
    
    
    %% bootstrap c values (average dist to centroid) and  make plot with confidence intervals
%     cis = cellfun(@(y) bootci(2000,@calcMeanDistToCent,y'),dataMatAllMeas,'uniformoutput',0);
 
%     [cs,values] = cellfun(@(y) calcMeanDistToCent(y'),dataMatAllMeas,'uniformoutput',0);
%     setAxis('on',0.75);
%     hold on
    
    %arrayfun(@(x) scatter(x,cs(x),50,toPlot.info.color{x},'filled'),1:length(cs));
%     arrayfun(@(x) scatter(x,values,50,toPlot.info.color{x},'filled'),1:length(cs)); 
%     arrayfun(@(x) errorbar(x,cs(x),cis{x}(1),cis{x}(2),'color',toPlot.info.color{x}),1:length(cs));
%     ylabel('Mean Euclidean Distance to Centroid Per Group')
%     labels = toPlot.info.names;
    
%     set(gca,'XTick',1:length(cs));
%     set(gca,'XTickLabel',labels,'FontSize',10);
%     
%     %scatter(1:length(cs),cs);
%     % errorbar(cs,cis);
%     saveas(gcf,[ip.Results.OutputDirectory filesep 'DistanceToCentroidPerGroup.fig']);
%     saveas(gcf,[ip.Results.OutputDirectory filesep 'DistanceToCentroidPerGroup.png']);
    
    %% Plot the distance to centroid control per group 
    
    
    %%
%     names = toPlot.info.names';
%     forCell  = num2cell([DB' Dunn' c1' c2']);
%     values = [names forCell];
%     % cell2dataset(values);
%     discrimValues =cell2table(values);
%     save([distanceDir filesep 'Table_Dunn_DB_C1_C2'],'discrimValues');
    
    if ip.Results.PCA
        
        % By PCA
        
        scoresPC = arrayfun(@(x) scores(grouping==x,1:forLine),1:nGroups,'uniformoutput',0);
        [DBPCA,DunnPCA] = arrayfun(@(x) whDiscriminationMeasures(scoresPC{x}',scoresPC{1}'),1:nCond);
        
        
        forCellPCA  = num2cell([DBPCA' DunnPCA']);
        valuesPCA = [names forCellPCA];
        % cell2dataset(values);
        discrimValuesPCA =cell2table(valuesPCA);
        save([distanceDir filesep 'Table_Dunn_DB_PCA'],'discrimValuesPCA','scoresPC');
    end
        
end


end
function [meanDistToCent,genesDist] = calcMeanDistToCent(featsVect)
%      featsVect = rxc double array 
%      where r is the number of features and c is the number of
%      observations (often cells/movies) of the perturbation condition
 % get the mean for each feature over all observations 
  meanGene = nanmean(featsVect,2)'; % column is 
  genesDist = pdist2(featsVect',meanGene); % get the eucledian distance between each observation and the centroid 
  meanDistToCent = nanmean(genesDist); 
end 

function [medianDistToCent,genesDist] = calcMedianDistToCent(featsVect)
%      featsVect = rxc double array 
%      where r is the number of features and c is the number of
%      observations (often cells/movies) of the perturbation condition
 % get the mean for each feature over all observations 
  medianGene = nanmean(featsVect,2)'; % column is 
  genesDist = pdist2(featsVect',medianGene); % get the eucledian distance between each observation and the centroid 
  medianDistToCent = nanmedian(genesDist); 
end 



%% 
