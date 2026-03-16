classdef FocalAdhesionPipeline < handle
    % FocalAdhesionPipeline - Wrapper for focal adhesion analysis pipeline
    % Unique name to avoid conflict with u-track's FocalAdhesionPackage
    %
    % Usage:
    %   pkg = FocalAdhesionPipeline('/path/to/Combined_Movie.tif');
    %   pkg.run();
    %
    % Options:
    %   pkg = FocalAdhesionPipeline(path, 'optimizationMode', 'cpu');
    %   pkg = FocalAdhesionPipeline(path, 'optimizationMode', 'gpu');
    %   pkg = FocalAdhesionPipeline(path, 'optimizationMode', 'none');
    %
    
    properties
        tiffPath
        MD
        FAPackage
        
        % Parameters (same defaults as base code)
        pixelSize = 108
        timeInterval = 10
        numAperture = 1.49
        camBitdepth = 16
        fluorophore = 'rfp'
        
        % Options
        forceRerun = true
        
        % Optimization modes for Step 7 (track processing):
        %   'none'     - Original code (~4 hours)
        %   'cpu'      - CPU-parallelized, exact output (~5-10 min)
        %   'gpu'      - GPU moment-based, fastest but different output (~1-2 min)
        optimizationMode = 'cpu'
        
        % Legacy support (maps to optimizationMode)
        parallelOptimization = true  % true='cpu', false='none'
        
        % Classification parameters
        startingDistG1G2 = 4  % Distance from cell edge for G1/G2 classification (µm)
        
        % Step 11 options
        runStep11 = true  % Enable/disable Step 11 (Initial Rise Time Lag Calculation)
        
        % Results
        metrics struct
        success = false
        timingInfo struct  % Store timing for each step
    end
    
    methods
        function obj = FocalAdhesionPipeline(tiffPath, varargin)
            obj.tiffPath = tiffPath;
            obj.timingInfo = struct();
            
            % Parse optional parameters
            for i = 1:2:length(varargin)
                if isprop(obj, varargin{i})
                    obj.(varargin{i}) = varargin{i+1};
                % Support legacy 'gpuOptimization' parameter name
                elseif strcmpi(varargin{i}, 'gpuOptimization')
                    if varargin{i+1}
                        obj.optimizationMode = 'cpu';
                    else
                        obj.optimizationMode = 'none';
                    end
                    obj.parallelOptimization = varargin{i+1};
                end
            end
            
            % Sync parallelOptimization with optimizationMode
            if obj.parallelOptimization && strcmp(obj.optimizationMode, 'none')
                obj.optimizationMode = 'cpu';
            elseif ~obj.parallelOptimization && ~strcmp(obj.optimizationMode, 'none')
                obj.optimizationMode = 'none';
            end
        end
        
        function run(obj)
            totalStartTime = tic;
            
           %% STEP 0: SETUP MOVIEDATA (exact copy from base code)
            stepStart = tic;
            channel = Channel(obj.tiffPath);   % obj.tiffPath is now a folder path
            analysisDir = fullfile(obj.tiffPath, 'analysis');
            if ~exist(analysisDir, 'dir'), mkdir(analysisDir); end
            MD = MovieData(channel, analysisDir);
            MD.outputDirectory_ = analysisDir;
            MD.movieDataPath_ = analysisDir;
            MD.movieDataFileName_ = 'movieData.mat';
             
            MD.pixelSize_ = obj.pixelSize;
            MD.timeInterval_ = obj.timeInterval;
            MD.numAperture_ = obj.numAperture;
            MD.camBitdepth_ = obj.camBitdepth;
            
            fluorProps = getFluorPropStruct();
            idx = find(strcmpi({fluorProps.name}, obj.fluorophore));
            if ~isempty(idx)
                MD.channels_(1).fluorophore_ = fluorProps(idx).name;
                MD.channels_(1).emissionWavelength_ = fluorProps(idx).lambda_em * 1e9;
            end
            
            MD.sanityCheck();
            MD.save();
            obj.MD = MD;
            
            disp('=== MovieData Setup Complete ===');
            disp(['  Folder: ' obj.tiffPath]);
            disp(['  Frames: ' num2str(MD.nFrames_)]);
            disp(['  Image size: ' num2str(MD.imSize_(1)) 'x' num2str(MD.imSize_(2))]);
            disp(['  Optimization mode: ' obj.optimizationMode]);
            disp(['  G1/G2 edge distance: ' num2str(obj.startingDistG1G2) ' µm']);
            obj.timingInfo.setup = toc(stepStart);
            %% STEP 0.5: ADD FOCAL ADHESION PACKAGE
            iPack = MD.getPackageIndex('FocalAdhesionPackage', 1, false);
            if isempty(iPack)
                MD.addPackage(FocalAdhesionPackage(MD));
                iPack = MD.getPackageIndex('FocalAdhesionPackage');
            end
            FAPackage = MD.getPackage(iPack);
            obj.FAPackage = FAPackage;
            
            %% STEP 1: STAGE DRIFT CORRECTION
            stepStart = tic;
            disp('=== Step 1: Stage Drift Correction ===');
            if isempty(FAPackage.processes_{1})
                FAPackage.createDefaultProcess(1);
            end
            proc1 = FAPackage.processes_{1};
            params1 = proc1.funParams_;
            params1.ChannelIndex = 1;
            params1.referenceFrameNum = 1;
            params1.usfac = 20;
            proc1.setPara(params1);
            MD.save();
            proc1.run();
            MD.save();
            obj.timingInfo.step1_driftCorrection = toc(stepStart);
            fprintf('  Step 1 completed in %.1f sec\n', obj.timingInfo.step1_driftCorrection);
            
            %% STEP 2: THRESHOLDING
            stepStart = tic;
            disp('=== Step 2: Thresholding ===');
            if isempty(FAPackage.processes_{2})
                FAPackage.createDefaultProcess(2);
            end
            proc2 = FAPackage.processes_{2};
            params2 = proc2.funParams_;
            params2.ChannelIndex = 1;
            params2.MethodIndx = 1;              % 1=MinMax, 2=Otsu, 3=Rosin, 4=Gradient-based
            params2.GaussFilterSigma = 2;        % Gaussian filter with 2 pixel std deviation
            proc2.setPara(params2);
            MD.save();
            proc2.run();
            MD.save();
            obj.timingInfo.step2_thresholding = toc(stepStart);
            fprintf('  Step 2 completed in %.1f sec\n', obj.timingInfo.step2_thresholding);

            
                        
            %% STEP 3: MASK REFINEMENT
            stepStart = tic;
            disp('=== Step 3: Mask Refinement ===');
            if isempty(FAPackage.processes_{3})
                FAPackage.createDefaultProcess(3);
            end
            proc3 = FAPackage.processes_{3};
            params3 = proc3.funParams_;
            params3.ChannelIndex = 1;
            params3.MaskCleanUp = 1;
            params3.MinimumSize = 10;
            params3.ClosureRadius = 3;
            params3.OpeningRadius = 0;
            params3.ObjectNumber = 1;
            params3.FillHoles = 1;
            proc3.setPara(params3);
            MD.save();
            proc3.run();
            MD.save();
            obj.timingInfo.step3_maskRefinement = toc(stepStart);
            fprintf('  Step 3 completed in %.1f sec\n', obj.timingInfo.step3_maskRefinement);
            
            %% STEP 4: POINT SOURCE DETECTION
            stepStart = tic;
            disp('=== Step 4: Point Source Detection ===');
            if isempty(FAPackage.processes_{4})
                FAPackage.createDefaultProcess(4);
            end
            proc4 = FAPackage.processes_{4};
            params4 = proc4.funParams_;
            params4.ChannelIndex = 1;
            params4.MaskChannelIndex = 1;
            params4.MaskProcessIndex = 3;
            params4.alpha = 0.05;
            params4.maskRadius = 40;
            params4.Mode = {'xyAc'};
            params4.FitMixtures = 0;
            params4.MaxMixtures = 5;
            params4.RedundancyRadius = 0.25;
            params4.UseIntersection = 1;
            proc4.setPara(params4);
            MD.save();
            proc4.run();
            MD.save();
            obj.timingInfo.step4_detection = toc(stepStart);
            fprintf('  Step 4 completed in %.1f sec\n', obj.timingInfo.step4_detection);
            
            %% STEP 5: TRACKING
            stepStart = tic;
            disp('=== Step 5: Tracking ===');
            if isempty(FAPackage.processes_{5})
                FAPackage.createDefaultProcess(5);
            end
            proc5 = FAPackage.processes_{5};
            params5 = proc5.funParams_;
            params5.ChannelIndex = 1;
            params5.probDim = 2;
            params5.verbose = 1;
            params5.gapCloseParam.timeWindow = 3;
            params5.gapCloseParam.minTrackLen = 1;
            params5.gapCloseParam.diagnostics = 0;
            params5.costMatrices(1).parameters.linearMotion = 2;
            params5.costMatrices(1).parameters.minSearchRadius = 1;
            params5.costMatrices(1).parameters.maxSearchRadius = 4;
            params5.costMatrices(1).parameters.brownStdMult = 3;
            params5.costMatrices(1).parameters.useLocalDensity = 1;
            params5.costMatrices(1).parameters.nnWindow = 3;
            params5.costMatrices(2).parameters.linearMotion = 2;
            params5.costMatrices(2).parameters.minSearchRadius = 1;
            params5.costMatrices(2).parameters.maxSearchRadius = 4;
            params5.costMatrices(2).parameters.maxAngleVV = 45;
            params5.costMatrices(2).parameters.gapPenalty = 1.5;
            params5.costMatrices(2).parameters.ampRatioLimit = [0.5 2];
            proc5.setPara(params5);
            MD.save();
            proc5.run();
            MD.save();
            obj.timingInfo.step5_tracking = toc(stepStart);
            fprintf('  Step 5 completed in %.1f sec\n', obj.timingInfo.step5_tracking);
            
            %% STEP 6: FA SEGMENTATION
            stepStart = tic;
            disp('=== Step 6: FA Segmentation ===');
            if isempty(FAPackage.processes_{6})
                FAPackage.createDefaultProcess(6);
            end
            proc6 = FAPackage.processes_{6};
            params6 = proc6.funParams_;
            params6.ChannelIndex = 1;
            params6.SteerableFilterSigma = 200;
            params6.OpeningRadiusXY = 0;
            params6.OpeningHeightT = 10;
            params6.MinVolTime = 1;
            proc6.setPara(params6);
            MD.save();
            proc6.run();
            MD.save();
            obj.timingInfo.step6_segmentation = toc(stepStart);
            fprintf('  Step 6 completed in %.1f sec\n', obj.timingInfo.step6_segmentation);
            
            %% STEP 7: ADHESION ANALYSIS
            stepStart = tic;
            disp('=== Step 7: Adhesion Analysis ===');
            if isempty(FAPackage.processes_{7})
                FAPackage.createDefaultProcess(7);
            end
            proc7 = FAPackage.processes_{7};
            
            if obj.forceRerun
                outputDir7 = proc7.funParams_.OutputDirectory;
                if exist(outputDir7, 'dir')
                    rmdir(outputDir7, 's');
                end
            end
            
            params7 = proc7.funParams_;
            params7.ChannelIndex = 1;
            params7.SegCellMaskProc = 3;
            params7.detectedNAProc = 4;
            params7.trackFAProc = 5;
            params7.FAsegProc = 6;
            params7.ApplyCellSegMask = 1;
            params7.reTrack = 1;
            
            % Set optimization method (new parameter name)
            % 'none' = original code (readIntensityFromTracks)
            % 'cpu'  = CPU-parallelized (readIntensityFromTracks_cpuParallel)
            % 'gpu'  = GPU moment-based (readIntensityFromTracks_gpuMomentBased)
            params7.optimization_method = obj.optimizationMode;
            switch obj.optimizationMode
                case 'none'
                    fprintf('  Using: Original code (slow)\n');
                case 'cpu'
                    fprintf('  Using: CPU-parallelized (readIntensityFromTracks_cpuParallel)\n');
                case 'gpu'
                    fprintf('  Using: GPU moment-based (readIntensityFromTracks_gpuMomentBased)\n');
                otherwise
                    params7.optimization_method = 'cpu';
                    fprintf('  Using: CPU-parallelized (default)\n');
            end
            
            params7.matchWithFA = 1;
            params7.onlyEdge = 0;
            params7.getEdgeRelatedFeatures = 1;
            params7.minLifetime = 3;
            params7.bandwidthNA = 4;
            params7.minFALengthMicron = 2;
            params7.backupOldResults = 0;
            params7.showAllTracks = 0;
            params7.plotEachTrack = 0;
            proc7.setPara(params7);
            MD.save();
            proc7.run();
            MD.save();
            obj.timingInfo.step7_analysis = toc(stepStart);
            fprintf('  Step 7 completed in %.1f sec\n', obj.timingInfo.step7_analysis);
            
            %% STEP 8: ADHESION CLASSIFICATION
            stepStart = tic;
            disp('=== Step 8: Adhesion Classification ===');
            if isempty(FAPackage.processes_{8})
                FAPackage.createDefaultProcess(8);
            end
            proc8 = FAPackage.processes_{8};
            
            if obj.forceRerun
                outputDir8 = proc8.funParams_.OutputDirectory;
                if exist(outputDir8, 'dir')
                    rmdir(outputDir8, 's');
                end
            end
            
            params8 = proc8.funParams_;
            params8.ChannelIndex = 1;
            params8.SegCellMaskProc = 3;
            params8.detectedNAProc = 4;
            params8.trackFAProc = 5;
            params8.FAsegProc = 6;
            params8.ApplyCellSegMask = 1;
            params8.useAutomaticallySelectedData = 1;
            params8.manualLabeling = 0;
            params8.useSimpleClassification = 0;
            params8.startingDist = obj.startingDistG1G2;
            params8.backupOldResults = 0;
            proc8.setPara(params8);
            MD.save();
            
            % Verify parameter was set correctly
            verifyParams = proc8.funParams_;
            fprintf('  [Verify] startingDist set to: %d µm\n', verifyParams.startingDist);
            
            proc8.run();
            MD.save();
            
            % Verify parameter after run (in case it was overwritten)
            finalParams = proc8.funParams_;
            if finalParams.startingDist ~= obj.startingDistG1G2
                warning('startingDist changed during run! Expected %d, got %d', obj.startingDistG1G2, finalParams.startingDist);
            end
            obj.timingInfo.step8_classification = toc(stepStart);
            fprintf('  Step 8 completed in %.1f sec\n', obj.timingInfo.step8_classification);
            
            %% STEP 11: INITIAL RISE TIME LAG CALCULATION (Optional)
            if obj.runStep11
                stepStart = tic;
                disp('=== Step 11: Initial Rise Time Lag Calculation ===');
                
                try
                    % Get or create Step 11 process
                    if isempty(FAPackage.processes_{11})
                        FAPackage.createDefaultProcess(11);
                    end
                    proc11 = FAPackage.processes_{11};
                    
                    % Run Step 11 - this will save MAT files to data/ subdirectory
                    calculateInitialRiseTimeLagFromTracks(MD);
                    MD.save();
                    
                    % Verify critical output files exist
                    step11OutputDir = proc11.funParams_.OutputDirectory;
                    step11DataDir = fullfile(step11OutputDir, 'data');
                    
                    requiredFiles = {'assemRate.mat', 'disassemRate.mat', 'idGroups.mat'};
                    missingFiles = {};
                    for f = 1:length(requiredFiles)
                        if ~exist(fullfile(step11DataDir, requiredFiles{f}), 'file')
                            missingFiles{end+1} = requiredFiles{f};
                        end
                    end
                    
                    if ~isempty(missingFiles)
                        error('Step 11 failed: Missing output files: %s', strjoin(missingFiles, ', '));
                    end
                    
                    obj.timingInfo.step11_timeLag = toc(stepStart);
                    fprintf('  Step 11 completed in %.1f sec\n', obj.timingInfo.step11_timeLag);
                    fprintf('  Output directory: %s\n', step11DataDir);
                    
                catch ME
                    error('Step 11 (Initial Rise Time Lag Calculation) failed: %s', ME.message);
                end
            else
                disp('=== Step 11: Skipped (runStep11 = false) ===');
            end
            
            obj.success = true;
            obj.timingInfo.total = toc(totalStartTime);
            
            disp('=== Pipeline Complete ===');
            fprintf('  Total time: %.1f sec (%.1f min)\n', obj.timingInfo.total, obj.timingInfo.total/60);
            obj.printTimingSummary();
        end
        
        function printTimingSummary(obj)
            % Print timing breakdown
            disp('--- Timing Summary ---');
            if isfield(obj.timingInfo, 'step1_driftCorrection')
                fprintf('  Step 1 (Drift Correction):  %6.1f sec\n', obj.timingInfo.step1_driftCorrection);
            end
            if isfield(obj.timingInfo, 'step2_thresholding')
                fprintf('  Step 2 (Thresholding):      %6.1f sec\n', obj.timingInfo.step2_thresholding);
            end
            if isfield(obj.timingInfo, 'step3_maskRefinement')
                fprintf('  Step 3 (Mask Refinement):   %6.1f sec\n', obj.timingInfo.step3_maskRefinement);
            end
            if isfield(obj.timingInfo, 'step4_detection')
                fprintf('  Step 4 (Detection):         %6.1f sec\n', obj.timingInfo.step4_detection);
            end
            if isfield(obj.timingInfo, 'step5_tracking')
                fprintf('  Step 5 (Tracking):          %6.1f sec\n', obj.timingInfo.step5_tracking);
            end
            if isfield(obj.timingInfo, 'step6_segmentation')
                fprintf('  Step 6 (Segmentation):      %6.1f sec\n', obj.timingInfo.step6_segmentation);
            end
            if isfield(obj.timingInfo, 'step7_analysis')
                fprintf('  Step 7 (Analysis):          %6.1f sec  <-- Optimized\n', obj.timingInfo.step7_analysis);
            end
            if isfield(obj.timingInfo, 'step8_classification')
                fprintf('  Step 8 (Classification):    %6.1f sec\n', obj.timingInfo.step8_classification);
            end
            if isfield(obj.timingInfo, 'step11_timeLag')
                fprintf('  Step 11 (Time Lag Calc):    %6.1f sec\n', obj.timingInfo.step11_timeLag);
            end
            disp('----------------------');
        end
        
        function metrics = getMetrics(obj)
            % Extract metrics after pipeline completes
            % Loads data from Step 11 MAT files (assembly/disassembly rates)
            % and from Step 7/8 for lifetimes
            
            metrics = struct();
            metrics.assembly = struct('G1', [], 'G4', [], 'G8', []);
            metrics.disassembly = struct('G1', [], 'G4', [], 'G8', []);
            metrics.lifetime = struct('G1', [], 'G2', [], 'totalFA', []);
            metrics.groupCounts = zeros(9, 1);
            
            if ~obj.success
                error('Pipeline has not completed successfully. Run pipeline first.');
            end
            
            try
                %% LOAD FROM STEP 11 MAT FILES (assembly/disassembly rates)
                if obj.runStep11
                    disp('Loading metrics from Step 11 MAT files...');
                    
                    % Locate Step 11 data directory
                    proc11 = obj.FAPackage.processes_{11};
                    step11OutputDir = proc11.funParams_.OutputDirectory;
                    step11DataDir = fullfile(step11OutputDir, 'data');
                    
                    % Load assembly rates (cell array: {G1, G2, ..., G9})
                    assemFile = fullfile(step11DataDir, 'assemRate.mat');
                    if exist(assemFile, 'file')
                        assemData = load(assemFile);
                        metrics.assembly.G1 = assemData.assemRate{1}(:);  % G1: NA Turnover
                        metrics.assembly.G4 = assemData.assemRate{4}(:);  % G4: Stalling
                        metrics.assembly.G8 = assemData.assemRate{8}(:);  % G8: Stable FA
                    else
                        error('Step 11 output not found: %s', assemFile);
                    end
                    
                    % Load disassembly rates
                    disassemFile = fullfile(step11DataDir, 'disassemRate.mat');
                    if exist(disassemFile, 'file')
                        disassemData = load(disassemFile);
                        metrics.disassembly.G1 = disassemData.disassemRate{1}(:);
                        metrics.disassembly.G4 = disassemData.disassemRate{4}(:);
                        metrics.disassembly.G8 = disassemData.disassemRate{8}(:);
                    else
                        error('Step 11 output not found: %s', disassemFile);
                    end
                    
                    disp('  Assembly/disassembly rates loaded from Step 11');
                else
                    warning('Step 11 was not run. Assembly/disassembly rates will be empty.');
                end
                
                %% LOAD LIFETIMES FROM STEP 7/8 (track files)
                disp('Loading lifetimes from track files...');
                
                proc7 = obj.FAPackage.processes_{7};
                proc8 = obj.FAPackage.processes_{8};
                
                analysisDir = proc7.funParams_.OutputDirectory;
                classDir = proc8.funParams_.OutputDirectory;
                
                % Load track folder and get format string
                metaFiles = dir(fullfile(analysisDir, '*metaTrackData.mat'));
                if isempty(metaFiles)
                    error('No metaTrackData.mat found in %s', analysisDir);
                end
                meta = load(fullfile(analysisDir, metaFiles(1).name));
                trackFolder = meta.metaTrackData.trackFolderPath;
                numTracks = meta.metaTrackData.numTracks;
                
                % Dynamic format string
                numDigits = floor(log10(numTracks)) + 1;
                fString = ['%0' num2str(numDigits) 'd'];
                
                % Load classification
                idsFile = dir(fullfile(classDir, '*idsClassified*.mat'));
                if isempty(idsFile)
                    error('No idsClassified file found in %s', classDir);
                end
                idsData = load(fullfile(classDir, idsFile(1).name));
                
                % Get group indices
                groupIndices = cell(9, 1);
                for g = 1:9
                    fieldName = ['idGroup' num2str(g)];
                    if isfield(idsData, fieldName)
                        val = idsData.(fieldName);
                        if islogical(val)
                            groupIndices{g} = find(val);
                        else
                            groupIndices{g} = val(:);
                        end
                        metrics.groupCounts(g) = length(groupIndices{g});
                    end
                end
                
                % Extract lifetimes per group
                timeInterval = obj.timeInterval;
                
                % G1: NA Turnover lifetimes
                for idx = groupIndices{1}'
                    track = obj.loadTrack(trackFolder, idx, fString);
                    if isempty(track), continue; end
                    metrics.lifetime.G1(end+1) = obj.getLifetime(track) * timeInterval;
                end
                
                % G2: Maturing lifetimes
                for idx = groupIndices{2}'
                    track = obj.loadTrack(trackFolder, idx, fString);
                    if isempty(track), continue; end
                    metrics.lifetime.G2(end+1) = obj.getLifetime(track) * timeInterval;
                end
                
                % G8 + G9: Total FA lifetimes
                for idx = [groupIndices{8}; groupIndices{9}]'
                    track = obj.loadTrack(trackFolder, idx, fString);
                    if isempty(track), continue; end
                    metrics.lifetime.totalFA(end+1) = obj.getLifetime(track) * timeInterval;
                end
                
                %% PRINT SUMMARY
                disp('--- Group Counts ---');
                groupNames = {'G1: NA Turnover', 'G2: Maturing', 'G3: Fading', ...
                              'G4: Stalling', 'G5: Uncertain', 'G6: Noisy', ...
                              'G7: Dissolving FC', 'G8: Strong Stable FA', 'G9: Weak Stable FA'};
                for g = 1:9
                    fprintf('  %s: %d\n', groupNames{g}, metrics.groupCounts(g));
                end
                
                disp('--- Metrics Summary ---');
                if ~isempty(metrics.assembly.G1)
                    fprintf('  G1 Assembly rate: %.3f ± %.3f (n=%d)\n', ...
                        mean(metrics.assembly.G1), std(metrics.assembly.G1), length(metrics.assembly.G1));
                end
                if ~isempty(metrics.lifetime.G1)
                    fprintf('  G1 Lifetime: %.1f ± %.1f sec (n=%d)\n', ...
                        mean(metrics.lifetime.G1), std(metrics.lifetime.G1), length(metrics.lifetime.G1));
                end
                if ~isempty(metrics.lifetime.totalFA)
                    fprintf('  FA Lifetime: %.1f ± %.1f sec (n=%d)\n', ...
                        mean(metrics.lifetime.totalFA), std(metrics.lifetime.totalFA), length(metrics.lifetime.totalFA));
                end
                
            catch ME
                error('Error extracting metrics: %s\n%s', ME.message, getReport(ME));
            end
            
            obj.metrics = metrics;
        end
    end
    
    methods (Access = private)
        function track = loadTrack(~, trackFolder, idx, fString)
            % FIXED: Use dynamic format string
            track = [];
            trackFile = fullfile(trackFolder, sprintf(['track' fString '.mat'], idx));
            if exist(trackFile, 'file')
                tmp = load(trackFile);
                track = tmp.curTrack;
            else
                % Try to find the file with glob pattern as fallback
                pattern = fullfile(trackFolder, sprintf('track*%d.mat', idx));
                files = dir(pattern);
                if ~isempty(files)
                    tmp = load(fullfile(trackFolder, files(1).name));
                    track = tmp.curTrack;
                end
            end
        end
        
        function lt = getLifetime(~, track)
            lt = 0;
            if isfield(track, 'presence')
                lt = sum(track.presence);
            elseif isfield(track, 'lifeTime')
                lt = track.lifeTime;
            elseif isfield(track, 'endingFrameExtra') && isfield(track, 'startingFrameExtra')
                lt = track.endingFrameExtra - track.startingFrameExtra + 1;
            end
        end
    end
end