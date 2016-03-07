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
ip.addParameter('plotByGroup',true); % will make 3-D scatter plots by group
%ip.addParameter('perFrame',true); Not implemented here: see other versions

% D Reduction Visualizations
ip.addParameter('PCA',true);
ip.addParameter('MDS',false);
%ip.addParameter('CriterionMDS','stress'); % default is nonmetric

ip.addParameter('cluster',true);

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


nGroups = numel(toPlot.info.names);
grouping = toPlot.info.grouping;

%% User can define which measurements(variables) to include if interactive turned on
if ip.Results.interactive
    paramSelect  = listSelectGUI(varNames,[],'move');
    varNames = varNames(paramSelect);
end
nVars = numel(varNames);
varNamesString = strrep(varNames, '_' ,'');
%% plot by frame : note this is going to need to be a different structure
% because lose information regarding the average per frame from the origina
%
% data format you want is frameX,cellX.
% if ip.Results.perFrame
%     for iGroup = 1:numel(toPlot.info.names)
%         dataMatC =  toPlot.(varNames{x}).dataMat{iGroup};
%     end
% else

%% Collect and normalize the measurement data : Default is to take a median value for each cell

dataDir = [ip.Results.OutputDirectory filesep 'DataFiles'];
if ~isdir(dataDir)
    mkdir(dataDir);
end

for iGroup  = 1:numel(toPlot.info.names)
    % transform data so that each row is a neurite and each column is a
    % set of measurements- this will be a per cell variable.
    perNeuriteStat = str2func(ip.Results.perNeuriteStatistic);
    dataMatC = arrayfun(@(x) perNeuriteStat(toPlot.(varNames{x}).dataMat{iGroup},1)',1:nVars,'uniformoutput',0);
    % put all the measurement values together (median measurement value for each neurite per column)
    % r indicates the neurite number;
    dataMatC = horzcat(dataMatC{:});
    
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

% Save the raw and normalized collected data in a csv.
% NOTE: need to change to table eventually to make compatible with later
% matlab versions!
dataPreNorm = mat2dataset(dataMatAllGroupsMeas,'ObsNames',obsNames,'varNames',varNames);
save([dataDir filesep 'OriginalValuesMedianPerNeuriteMovie'],'dataPreNorm');
export(dataPreNorm,'file',[dataDir filesep 'OriginalValuesMedianPerNeuriteMovie.csv']);
%writetable(tDataRaw,[ip.Results.OutputDirectory filesep 'OriginalValuesMedianPerNeuriteMovie.csv'],dataOriginal);

dataZScores = mat2dataset(dataFinal,'ObsNames',obsNames,'varNames',varNames);
export(dataZScores,'file',[dataDir filesep 'dataZScores.csv']);
save([dataDir filesep 'dataZScores'],'dataZScores');

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





if ip.Results.cluster
    clusterDir = [ip.Results.OutputDirectory filesep 'Cluster'];
    if ~isdir(clusterDir)
        mkdir(clusterDir);
    end
    
    colorClust{1} = [0 0 1]; % low blue
    colorClust{2} = [1 0 0]; % high red
     
    % Perform k-means clustering (All Values)
    [idx,cCenters]  = kmeans(plotValues,2,'replicates',20);
    indexMin = find(cCenters == min(cCenters));
    indexMax = find(cCenters == max(cCenters));
    %  % make it such that 2 is always high and 1 is always the low cluster
    sortedIdx = zeros(length(idx),1);
    sortedIdx(idx==indexMin) = 1;
    sortedIdx(idx==indexMax) = 2;
    groupingCluster = sortedIdx;
    
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
    criterion{2} = 'metricstress';
    criterion{3} = 'sammon';
    criterion{4} = 'strain';
    
    for iCrit = 1:numel(criterion)
        
        cOutMDS =  [outMDS filesep criterion{iCrit} ];
        if ~isdir(cOutMDS)
            mkdir(cOutMDS)
        end
        
        cOutMDSDia = [cOutMDS filesep 'Diagnostic'];
        
        if ~isdir(cOutMDSDia)
            mkdir(cOutMDSDia);
        end
        
        % perform the nonclassical multi-dimensional scaling
        [ y, stress, disparities] = mdscale(dissimilarities,2,'Criterion',criterion{iCrit});
       % if strcmpi(criterion{iCrit},'stress');
            save([cOutMDSDia filesep 'MDSResults' criterion{iCrit}],'y','stress','disparities');
        %else
         %   save([cOutMDSDia filesep 'MDSResults' criterion{iCrit}],'y','stress');
        %end
        axisLims = [floor(min(y(:,1))-1),ceil(max(y(:,1)))+1,floor(min(y(:,2))-1),ceil(max(y(:,2))+1)];
        
        %% Plot shepards plot
        setAxis('on');
        distances = pdist(y);
        [dum,ord] = sortrows([disparities(:) dissimilarities(:)]);
        plot(dissimilarities,distances,'bo', ...
            dissimilarities(ord),disparities(ord),'r.-', ...
            [0 25],[0 25],'k-')
        %transformed values of the original dissimilarities
        
        xlabel('Dissimilarities')
        ylabel('Distances/Disparities')
        legend({'Distances' 'Disparities' '1:1 Line'},...
            'Location','NorthWest');
        saveas(gcf,[cOutMDSDia filesep 'shepardsPlot.fig']);
        saveas(gcf,[cOutMDSDia filesep 'shepardsPlot.png']);
        
        close gcf
        
        
        %% Plot in MDS space By Group
        
        setAxis('on')
        
        arrayfun(@(x) scatter(y(grouping==x,1),...
            y(grouping==x,2),50,toPlot.info.color{x},'filled'),1:nGroups);
        xlabel('MDS1');
        ylabel('MDS2');
        saveas(gcf,[cOutMDS filesep 'MDS_2DScatterByGroupColor.fig']);
        saveas(gcf,[cOutMDS filesep 'MDS_2DScatterByGroupColor.eps'],'psc2');
        saveas(gcf,[cOutMDS filesep 'MDS_2DScatterByGroupColor.png']);
        
        close gcf

        %% plot control by cluster
        
        cOutMDSClust = [cOutMDS filesep 'cluster'];
        if ~isdir(cOutMDSClust)
            mkdir(cOutMDSClust)
        end
        
        controlValues = y(grouping==1,:);
        %      groupingCluster = toPlot.info.groupingIdxClust;
        groupingControl = groupingCluster(grouping==1);
       
        nClusts = length(unique(groupingCluster));
        setAxis('on')
        hold on
        arrayfun(@(x) scatter(controlValues(groupingControl==x,1),...
            controlValues(groupingControl==x,2),50, colorClust{x},'filled'),1:nClusts);
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
        name{2} = 'High';
        for iClust = 1:2
            setAxis('on')
            
            hold on
            
            %        scatter(controlValues(groupingControl==iClust,1),...
            %            controlValues(groupingControl==iClust,2),50, colorClust{iClust},'filled');
            
            arrayfun(@(x) scatter(controlValues(groupingControl==x,1),...
                controlValues(groupingControl==x,2),50, colorClust{x},'filled'),1:nClusts);
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
            y(groupingCluster==x,2),50, colorClust{x},'filled'),1:nClusts);
        
        
        %arrayfun(@(x) text(y(x,1),y(x,2),obsNames{x}),1:size(y(:,1)))
        xlabel('MDS1');
        ylabel('MDS2');
        
        axis(axisLims);
        
        saveas(gcf,[cOutMDSClust filesep 'MDS_2DScatterByGroupClusterAll.fig']);
        saveas(gcf,[cOutMDSClust filesep 'MDS_2DScatterByGroupClusterAll.eps'],'psc2');
        saveas(gcf,[cOutMDSClust filesep 'MDS_2DScatterByGroupClusterAll.png']);
        close gcf
        
        %% plot results by outgrowth scaleMap
        if ip.Results.plotByGroup
            %
            MDSGroupDir = [cOutMDS filesep 'MDSPerGroup'];
            if ~isdir(MDSGroupDir);
                mkdir(MDSGroupDir);
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
                %% plot results by group with picture 
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
                          
            %% plot for each group
        end % if ip.Results.plotByGroup
    end % if iCrit
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
    
    
    
    
    
    
    
    
    
    
    %% Discrimination Metrics
    nCond = numel(dataMatAllMeas);
    [DB,Dunn,SC,c1,c2] = arrayfun(@(x) whDiscriminationMeasures(dataMatAllMeas{x}',dataMatAllMeas{1}'),1:nCond);
    
    
    %% bootstrap c values (average dist to centroid) and  make plot with confidence intervals
%     cis = cellfun(@(y) bootci(2000,@calcMeanDistToCent,y'),dataMatAllMeas,'uniformoutput',0);
 
%     [cs,values] = cellfun(@(y) calcMeanDistToCent(y'),dataMatAllMeas,'uniformoutput',0);
%     setAxis('on',0.75);
%     hold on
    
    %arrayfun(@(x) scatter(x,cs(x),50,toPlot.info.color{x},'filled'),1:length(cs));
%     arrayfun(@(x) scatter(x,values,50,toPlot.info.color{x},'filled'),1:length(cs)); 
    arrayfun(@(x) errorbar(x,cs(x),cis{x}(1),cis{x}(2),'color',toPlot.info.color{x}),1:length(cs));
    ylabel('Mean Euclidean Distance to Centroid Per Group')
    labels = toPlot.info.names;
    
    set(gca,'XTick',1:length(cs));
    set(gca,'XTickLabel',labels,'FontSize',10);
    
    %scatter(1:length(cs),cs);
    % errorbar(cs,cis);
    saveas(gcf,[ip.Results.OutputDirectory filesep 'DistanceToCentroidPerGroup.fig']);
    saveas(gcf,[ip.Results.OutputDirectory filesep 'DistanceToCentroidPerGroup.png']);
    
    %% Plot the distance to centroid control per group 
    
    
    %%
    names = toPlot.info.names';
    forCell  = num2cell([DB' Dunn' c1' c2']);
    values = [names forCell];
    % cell2dataset(values);
    discrimValues =cell2table(values);
    save([distanceDir filesep 'Table_Dunn_DB_C1_C2'],'discrimValues');
    
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

%% 
