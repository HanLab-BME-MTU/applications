function [ output_args ] = GCAAnalysisPartitionTimeSeriesMovie(MD,varargin)
%GCAAnalysisPartitionTimeSeriesMovie 
% Movie Wrapper for GCAAnalysisPartitionTimeSeries
% Run the pause finder for outgrowth trajectories 
% make a nice outgrowth overlay (require fixing of the spy function )
% collect local params. 

% PARAMETERS

% defaultOutDir = [MD.outputDirectory_ filesep 'MEASUREMENT_EXTRACTION' filesep ... 
%         'Partition_Outgrowth_Trajectory']; 
 defaultOutDir = [MD.outputDirectory_ filesep 'MEASUREMENT_EXTRACTION' filesep ... 
         'Partition_Outgrowth_Trajectory_WithGrowthTimes'];
 defaultMeasFolder =  ['SegmentationPackage/StepsToReconstructTestBugFix20160426/'... 
     'GCAMeasurementExtraction_test20160510/WholeNeurite']; % currently the piece after GrowthConeAnalyzer 
 % and before Descriptor...
     
%%Input check
ip = inputParser;
ip.addParameter('OutputDirectory',defaultOutDir,@(x) ischar(x));
ip.addParameter('splineParam',0.01, @(x) isscalar(x) || isempty(x) ); 
ip.addParameter('threshPause',0.5, @isscalar); 
ip.addParameter('MeasurementFolder',defaultMeasFolder); 

ip.parse(varargin{:});

outDir = ip.Results.OutputDirectory; 
splineParam = ip.Results.splineParam; 
threshPause = ip.Results.threshPause; 
%%
outgrowthFile = [MD.outputDirectory_ filesep 'MEASUREMENT_EXTRACTION' ... 
        filesep 'GlobalFunctional' filesep 'neurite_outgrowth_measurements'... 
        filesep 'neuriteLengthOutput.mat'];
% load the neurite Length info 
if exist(outgrowthFile)~= 0 ;
    load(outgrowthFile);
    
    if ~isdir(outDir);
        mkdir(outDir);
    end
    % Here we will simply partition based on pausing in the neurite
    % trajectory (when the velocity of outgrowth reaches below a user
    % selected threshold)
    
    title = gcaGetNeuriteID(MD.outputDirectory_);
    
    globalMeas =  GCAfindPausingInNeuriteOutgrowthTrajectory(neuriteLength,'outPath',outDir,... 
        'forTitle',title,'splineParam',splineParam,'threshPause',threshPause) ;
    % load
    % load([MD.outputDirectory filesep 'PARAMETER_EXTRACTION_201503015' filesep 'localParams.mat']);
    % go on to partition the rest of the data if user desires..
    %  grouping = globalMeas.grouping;
    
    %% for now just searchFiles
    parameterDir = [MD.outputDirectory_ filesep ip.Results.MeasurementFolder];
    % for now just search files - redesign so that the parameters in the
    % future are more cleverly named
    
    % search all descriptor parameters.
    localParamFiles = searchFiles('meas_',[],[parameterDir filesep 'Descriptor'],1);
    
    paramNames = cellfun(@(x) strrep(x,'meas_',''),localParamFiles(:,1),'uniformoutput',0);
    paramNames = cellfun(@(x) strrep(x,'.mat',''),paramNames,'uniformoutput',0);
    
    % for now (as of 20151214 do NOT include spatial autocorrelation)
    toKeep = ~strcmpi(paramNames,'maxACFLagSpatial'); 
    paramNames = paramNames(toKeep); 
    localParamFiles = localParamFiles(toKeep,:); 
    % for now just create a structure
    for iParam = 1:size(paramNames,1)
        load([localParamFiles{iParam,2} filesep localParamFiles{iParam,1}]);
       
        measCR = reformatDataCell(measC);
        localMeas.(paramNames{iParam}) = measCR;
    end
       
    %% Group all non-marco parameters based on outgrowth and get the predictor and response values
    [localMeas,predictors,responses] = GCAAnalysisPartitionTimeSeries(localMeas,globalMeas); %
    save([outDir filesep 'localMeas.mat'],'localMeas');
    %% Group all marco veil params (persTime of local veil protrusion and velocity veil prot)
        
    % load marcos analysis where he identifies viable protrusion/retraction events;
    %    load([MD.outputDirectory_ filesep 'EdgeVelocityQuantification' filesep 'EdgeMotion.mat']);
    protDir = [MD.outputDirectory_ filesep 'protrusion_samples_ConstantNumber_windSize_5ReInit61']; 
    
    if ~exist(protDir)
        protDir = [MD.outputDirectory_ filesep 'protrusion_samples_ConstantNumber_windSize_5_ReInit62']; 
    end 
        
    load([protDir filesep 'EdgeVelocityQuantification_CutMovie_1' ...
        filesep 'EdgeMotion.mat']);
    
    analysisResults1 = analysisResults;
    clear analysisResults
    load([protDir filesep 'EdgeVelocityQuantification_CutMovie_2' ...
        filesep 'EdgeMotion.mat']);
    analysisResults2 = analysisResults;
    clear analysisResults
    
    % now partition the veil params - this is a little more tricky as the
    % veil params intrinsically not measured per frame but the event may span
    % several frames choose only events that are completed within the given
    % window.
    [localMeasVeil,predictorsVeil,forMovie] = GCAAnalysisPartitionVeilParamsByOutgrowth(analysisResults1,analysisResults2,globalMeas);
    % for now save the veil and the filoBranch parameters in separate
    % structures
    save([outDir filesep 'localMeasVeil.mat'],'localMeasVeil');
    if ~isempty(forMovie)
        movieDir = [outDir filesep 'forMovie'];
        if ~isdir(movieDir)
            mkdir(movieDir)
        end
        save([movieDir filesep 'forMovie.mat'],'forMovie');
        
    end
        
    %     [localMeasVeil2,predictorsVeil2] = GCAAnalysisPartitionVeilParamsByOutgrowthNew(analysisResults1,globalMeas);
    
    % Create the dataset array of predictors and responses (response will be
    % the last
    
    paramNamesVeil = fieldnames(localMeasVeil);
    forMultReg = [predictors predictorsVeil responses];
    %    forMultReg = [predictors responses];
    forMultRegCell = num2cell(forMultReg);
    responseName = {'Neurite_Velocity'};
    %paramNamesAll = [paramNames' responseName];
    measNamesAll = [paramNames' paramNamesVeil' responseName];
    forMultRegCell = [measNamesAll;forMultRegCell];
    
    
    % save the lifetimes of each event etc for filtering 
    lifetimes =  cellfun(@(x) length(x),globalMeas.outgrowth.groupedFrames)'; 
    maxVel =  cellfun(@(x) max(x), globalMeas.outgrowth.groupedVelSmoothed)'; 
    stateGrp = globalMeas.outgrowth.stateGrp; 
    
%     filterInfo.stateGrp  = stateGrp; 
%     filterInfo.maxVel = maxVel; 
%     filterInfo.lifetime = lifetimes; 
filterInfo = [stateGrp maxVel lifetimes]; 
    
    multRegDir =   [outDir filesep 'ForMultRegression'];
    if ~isdir(multRegDir)
        mkdir(multRegDir)
    end
    %save([multRegDir filesep 'forMultRegCell.mat'],'forMultRegCell');
    % save as a .mat file with the names
    save([multRegDir filesep 'forMultReg'],'forMultReg' ,'measNamesAll','filterInfo'); %
    
    % save state information 
    
    dataSet = cell2dataset(forMultRegCell);
    
    % get a name for the .csv file
    [~,num] = upDirectory(MD.outputDirectory_,2,1);
    [~,date] = upDirectory(MD.outputDirectory_,3,1);
    [~,name] = upDirectory(MD.outputDirectory_,4,1);
    
    
    
    export(dataSet,'File',[multRegDir filesep name '_' date '_' num '.csv'],'Delimiter',',');
    %  save([multRegDir filesep 'predictors.mat'],'predictors');
    % save([multRegDir filesep 'responses.mat'],'responses');
    
   
    
    
else
    display('Please Run the Neurite Outgrowth Measurements to Continues');
end









end
    



