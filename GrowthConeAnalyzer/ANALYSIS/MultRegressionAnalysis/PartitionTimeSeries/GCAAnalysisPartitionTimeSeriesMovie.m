function [ output_args ] = GCAAnalysisPartitionTimeSeriesMovie(MD)
%GCAAnalysisPartitionTimeSeriesMovie 
% Movie Wrapper for GCAAnalysisPartitionTimeSeries
% Run the pause finder for outgrowth trajectories 
% make a nice outgrowth overlay (require fixing of the spy function )
% collect local params. 

outgrowthFile = [MD.outputDirectory_ filesep 'PARAMETER_EXTRACTION' ... 
        filesep 'GlobalFunctional' filesep 'neurite_outgrowth_measurements'... 
        filesep 'neuriteLengthOutput.mat'];
% load the neurite Length info 
if exist(outgrowthFile)~= 0 ;
    load(outgrowthFile); 
    outDir = [MD.outputDirectory_ filesep 'PARAMETER_EXTRACTION' filesep ... 
        'Partition_Outgrowth_Trajectory']; 
    if ~isdir(outDir); 
        mkdir(outDir); 
    end 
    % Here we will simply partition based on pausing in the neurite
    % trajectory (when the velocity of outgrowth reaches below a user 
    % selected threshold)
    globalParams =  GCAfindPausingInNeuriteOutgrowthTrajectory(neuriteLength,'outPath',outDir) ;
    % load 
   % load([MD.outputDirectory filesep 'PARAMETER_EXTRACTION_201503015' filesep 'localParams.mat']); 
    % go on to partition the rest of the data if user desires.. 
   %  grouping = globalParams.grouping;
  
  %% for now just searchFiles  
 parameterDir = [MD.outputDirectory_ filesep 'PARAMETER_EXTRACTION_20150305' ]; 
 % for now just search files - redesign so that the parameters in the
 % future are more cleverly named 
 
 % search all descriptor parameters. 
 localParamFiles = searchFiles('param_',[],[parameterDir filesep 'Descriptor'],1); 
  
 paramNames = cellfun(@(x) strrep(x,'param_',''),localParamFiles(:,1),'uniformoutput',0); 
 paramNames = cellfun(@(x) strrep(x,'.mat',''),paramNames,'uniformoutput',0); 
 
 % for now just create a structure 
 for iParam = 1:size(localParamFiles,1)
     load([localParamFiles{iParam,2} filesep localParamFiles{iParam,1}]); 
     paramCR = reformatDataCell(paramC); 
     localParams.(paramNames{iParam}) = paramCR;   
 end 
 
 
  
  %% Group all non-marco parameters based on outgrowth and get the predictor and response values  
   [localParms,predictors,responses] = GCAAnalysisPartitionTimeSeries(localParams,globalParams); % 
   save([outDir filesep 'localParams.mat'],'localParams'); 
  %% Group all marco veil params (persTime of local veil protrusion and velocity veil prot)  
  
  % load marcos analysis where he identifies viable protrusion/retraction events; 
   load([MD.outputDirectory_ filesep 'EdgeVelocityQuantification' filesep 'EdgeMotion.mat']); 
  
  
   
   groupedFrames = globalParams.outgrowth.groupedFrames; 
   
   % now partition the veil params - this is a little more tricky as the 
   % veil params intrinsically not measured per frame but the event may span 
   % several frames choose only events that are completed within the given 
   % window. 
   [localParamsVeil,predictorsVeil] = GCAAnalysisPartitionVeilParamsByOutgrowth(analysisResults,globalParams); 
   
    
   
   % Create the dataset array of predictors and responses (response will be
   % the last 
   
   paramNamesVeil = fieldnames(localParamsVeil); 
   forMultReg = [predictors predictorsVeil responses];
   forMultRegCell = num2cell(forMultReg); 
   responseName = {'Neurite_Velocity'}; 
   paramNamesAll = [paramNames' paramNamesVeil' responseName];
   forMultRegCell = [paramNamesAll;forMultRegCell]; 
   
   
   multRegDir =   [outDir filesep 'ForMultRegression'];
   if ~isdir(multRegDir) 
       mkdir(multRegDir)
   end 
   %save([multRegDir filesep 'forMultRegCell.mat'],'forMultRegCell'); 
   % save as a .mat file with the names 
   save([multRegDir filesep 'forMultReg'],'forMultReg' ,'paramNamesAll'); %
   
   
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
    



