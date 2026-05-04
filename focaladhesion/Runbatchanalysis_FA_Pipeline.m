%% Batch Focal Adhesion Analysis using FocalAdhesionPipeline
% Processes all Control and MEKi movies from Run 7 dataset
% Uses FocalAdhesionPipeline (NOT FocalAdhesionPackage) to avoid name conflict
% Fluorescence: TIRF mCherry (default) | Time interval: 1s

clear; clc; close all;

%% Configuration
% Run 7 dataset
basePath = '/mnt/nas/Collaborations/Michille/Run 8/run8-raw_downloaded/run8';
%% Timing
timeInterval = 10;  % seconds

%% Dynamically discover all crop folders and their *_combined.tif files
fprintf('\nSearching for *_combined.tif files in crop folders...\n');
allCropDirs = genpath(basePath);
allCropDirs = strsplit(allCropDirs, pathsep);
allCropDirs = allCropDirs(~cellfun('isempty', allCropDirs));

controlMovies = {};
mekiMovies = {};
for idx = 1:length(allCropDirs)
    [~, folderName] = fileparts(allCropDirs{idx});
    if ~contains(folderName, 'crop', 'IgnoreCase', true)
        continue;
    end
    
    % Check if folder contains any tif frames
    frameFiles = dir(fullfile(allCropDirs{idx}, '*.tif'));
    if isempty(frameFiles)
        continue;
    end
    
    % The folder itself is now the "movie"
    moviePath = allCropDirs{idx};
    
    % Classify condition based on path
    if contains(moviePath, 'control', 'IgnoreCase', true)
        controlMovies{end+1} = moviePath;
    elseif contains(moviePath, 'MEKi', 'IgnoreCase', true)
        mekiMovies{end+1} = moviePath;
    else
        fprintf('  WARNING: Could not classify: %s\n', moviePath);
    end
end
controlMovies = controlMovies(:);
mekiMovies = mekiMovies(:);

fprintf('Found %d Control movies and %d MEKi movies\n', length(controlMovies), length(mekiMovies));

%% TEST MODE - Process only 1 Control + 1 MEKi for quick feedback
testMode = false;  % Set to false for full batch processing

if testMode
    fprintf('========================================\n');
    fprintf('  TEST MODE ENABLED\n');
    fprintf('  Processing 1 Control + 1 MEKi only\n');
    fprintf('========================================\n\n');
    
    if ~isempty(controlMovies), controlMovies = controlMovies(1); end
    if ~isempty(mekiMovies), mekiMovies = mekiMovies(1); end
end

outputDir = fullfile(basePath, 'BatchAnalysisResults');
if ~exist(outputDir, 'dir'), mkdir(outputDir); end




%% Options
skipProcessed = false;  % Set to true to skip movies with existing results
forceRerun = true;      % Force rerun of adhesion analysis steps
suppressWarnings = true;  % Suppress MATLAB warnings during processing

% Optimization mode for track processing:
%   'none' - Original code (~4 hours per movie)
%   'cpu'  - CPU-parallelized with parfor (~10-15 min per movie) - EXACT OUTPUT
%   'gpu'  - GPU moment-based fitting (~2-5 min per movie) - DIFFERENT OUTPUT (scientifically valid)
optimizationMode = 'cpu';  % Recommended: 'cpu' for exact output, 'gpu' for maximum speed

%% Suppress warnings if requested
if suppressWarnings
    % Save current warning state
    originalWarningState = warning;
    
    % Turn off common warnings
    warning('off', 'all');  % Suppress all warnings
    
    % Or selectively suppress specific warnings:
    % warning('off', 'MATLAB:interp1:NaNstrip');
    % warning('off', 'MATLAB:nearlySingularMatrix');
    % warning('off', 'curvefit:fit:noStartPoint');
    % warning('off', 'stats:nlinfit:IllConditionedJacobian');
    % warning('off', 'MATLAB:singularMatrix');
    % warning('off', 'parallel:gpu:device:DeviceDeprecated');
    
    fprintf('  Warnings suppressed\n');
end

%% Initialize results storage
resultsFile = fullfile(outputDir, 'BatchResults.mat');

if skipProcessed && exist(resultsFile, 'file')
    load(resultsFile, 'results', 'processedMovies', 'movieTimings');
    fprintf('Loaded existing results with %d processed movies\n', length(processedMovies));
else
    results.control = initMetricsStruct();
    results.meki = initMetricsStruct();
    processedMovies = {};
    movieTimings = struct('movie', {}, 'condition', {}, 'time', {}, 'success', {});
end

%% Process all movies
allMovies = [controlMovies; mekiMovies];
allConditions = [repmat({'Control'}, length(controlMovies), 1); repmat({'MEKi'}, length(mekiMovies), 1)];

totalStartTime = tic;
successCount = 0;
failCount = 0;

fprintf('\n');
fprintf('========================================\n');
fprintf('  BATCH FOCAL ADHESION ANALYSIS\n');
fprintf('========================================\n');
fprintf('  Total movies: %d (%d Control, %d MEKi)\n', length(allMovies), length(controlMovies), length(mekiMovies));
fprintf('  Optimization mode: %s\n', optimizationMode);
fprintf('  Skip processed: %s\n', mat2str(skipProcessed));
fprintf('  Warnings suppressed: %s\n', mat2str(suppressWarnings));
fprintf('========================================\n\n');

for i = 1:length(allMovies)
    moviePath = allMovies{i};
    condition = allConditions{i};
    
    % Extract short name for display
    [~, parentDir] = fileparts(fileparts(moviePath));
    [~, grandParentDir] = fileparts(fileparts(fileparts(moviePath)));
    shortName = [grandParentDir '/' parentDir];
    
    fprintf('\n');
    fprintf('+--------------------------------------------------------------+\n');
    fprintf('¦ Movie %d/%d: %-20s [%s]\n', i, length(allMovies), shortName, condition);
    fprintf('+--------------------------------------------------------------+\n');
    
    % Check if already processed
    if skipProcessed && ismember(moviePath, processedMovies)
        fprintf('  [SKIPPED] Already processed\n');
        continue;
    end
    
    % Check if file exists
    if ~exist(moviePath, 'dir')
        fprintf('  [SKIPPED] File not found: %s\n', moviePath);
        failCount = failCount + 1;
        continue;
    end
    
    movieStartTime = tic;
    
    try
        % Run pipeline using FocalAdhesionPipeline with Step 11 enabled
        pkg = FocalAdhesionPipeline(moviePath, ...
            'forceRerun', forceRerun, ...
            'optimizationMode', optimizationMode, ...
            'timeInterval', timeInterval, ...
            'runStep11', true);  % Enable Step 11 for assembly/disassembly rate calculation
        pkg.run();
        
        % Get metrics (loads from Step 11 MAT files)
        metrics = pkg.getMetrics();
        
        % Append to results
        if strcmpi(condition, 'Control')
            results.control = appendMetrics(results.control, metrics);
        else
            results.meki = appendMetrics(results.meki, metrics);
        end
        
        movieTime = toc(movieStartTime);
        successCount = successCount + 1;
        
        % Record timing
        movieTimings(end+1) = struct('movie', moviePath, 'condition', condition, ...
                                      'time', movieTime, 'success', true);
        processedMovies{end+1} = moviePath;
        
        fprintf('\n');
        fprintf('  ? [SUCCESS] Completed in %.1f min\n', movieTime/60);
        
    catch ME
        movieTime = toc(movieStartTime);
        failCount = failCount + 1;
        
        % Record failure
        movieTimings(end+1) = struct('movie', moviePath, 'condition', condition, ...
                                      'time', movieTime, 'success', false);
        
        fprintf('\n');
        fprintf('  ? [FAILED] %s\n', ME.message);
        fprintf('    Stack trace:\n');
        for k = 1:min(3, length(ME.stack))
            fprintf('      %s (line %d)\n', ME.stack(k).name, ME.stack(k).line);
        end
        
        % Save error log
        errorLogFile = fullfile(outputDir, sprintf('error_movie%d_%s.txt', i, datestr(now, 'yyyymmdd_HHMMSS')));
        fid = fopen(errorLogFile, 'w');
        fprintf(fid, 'Movie: %s\n', moviePath);
        fprintf(fid, 'Error: %s\n\n', ME.message);
        fprintf(fid, '%s\n', getReport(ME, 'extended'));
        fclose(fid);
    end
    
    % Save intermediate results after each movie
    save(resultsFile, 'results', 'processedMovies', 'movieTimings');
    fprintf('  Intermediate results saved\n');
    
    % Print running summary
    elapsed = toc(totalStartTime);
    avgTime = elapsed / i;
    remaining = avgTime * (length(allMovies) - i);
    fprintf('  Progress: %d/%d | Success: %d | Failed: %d | ETA: %.1f min\n', ...
        i, length(allMovies), successCount, failCount, remaining/60);
end

%% Final Summary
totalTime = toc(totalStartTime);

fprintf('\n');
fprintf('========================================\n');
fprintf('  BATCH PROCESSING COMPLETE\n');
fprintf('========================================\n');
fprintf('  Total time: %.1f min (%.1f hours)\n', totalTime/60, totalTime/3600);
fprintf('  Success: %d/%d movies\n', successCount, length(allMovies));
fprintf('  Failed: %d movies\n', failCount);
fprintf('========================================\n');

% Print timing breakdown
if ~isempty(movieTimings)
    fprintf('\n--- Per-Movie Timing ---\n');
    successfulTimings = [movieTimings([movieTimings.success]).time];
    if ~isempty(successfulTimings)
        fprintf('  Average: %.1f min\n', mean(successfulTimings)/60);
        fprintf('  Min: %.1f min\n', min(successfulTimings)/60);
        fprintf('  Max: %.1f min\n', max(successfulTimings)/60);
    end
end

%% Generate plots
if successCount > 0
    fprintf('\nGenerating comparison plots...\n');
    try
        generatePlots(results, outputDir);
        fprintf('Plots saved to: %s\n', outputDir);
    catch ME
        fprintf('Error generating plots: %s\n', ME.message);
    end
end

%% Save combined data .mat file for later use
combinedDataFile = fullfile(outputDir, 'CombinedData_Run8.mat');

save(combinedDataFile, 'results', 'processedMovies', 'movieTimings', ...
     'controlMovies', 'mekiMovies', 'basePath', 'timeInterval', 'optimizationMode');
fprintf('\nCombined data saved to: %s\n', combinedDataFile);
fprintf('  Load later with: load(''%s'')\n', combinedDataFile);

% Print final data summary
printDataSummary(results);

%% Restore warnings
if suppressWarnings
    warning(originalWarningState);
    fprintf('\nWarnings restored to original state\n');
end

fprintf('\nDone!\n');

%% ========================================================================
%  HELPER FUNCTIONS
%  ========================================================================

function s = initMetricsStruct()
    s.assembly = struct('G1', [], 'G4', [], 'G8', [], 'All', []);
    s.disassembly = struct('G1', [], 'G4', [], 'G8', [], 'All', []);
    s.lifetime = struct('G1', [], 'G2', [], 'totalFA', [], 'All', []);
    s.groupCounts = [];  % Track group counts per movie
end

function r = appendMetrics(r, m)
    % Safely append metrics, handling empty/missing fields
    fields = {'G1', 'G4', 'G8'};
    for f = 1:length(fields)
        fn = fields{f};
        if isfield(m.assembly, fn) && ~isempty(m.assembly.(fn))
            r.assembly.(fn) = [r.assembly.(fn); m.assembly.(fn)(:)];
            r.assembly.All = [r.assembly.All; m.assembly.(fn)(:)];
        end
        if isfield(m.disassembly, fn) && ~isempty(m.disassembly.(fn))
            r.disassembly.(fn) = [r.disassembly.(fn); m.disassembly.(fn)(:)];
            r.disassembly.All = [r.disassembly.All; m.disassembly.(fn)(:)];
        end
    end
    
    lifeFields = {'G1', 'G2', 'totalFA'};
    for f = 1:length(lifeFields)
        fn = lifeFields{f};
        if isfield(m.lifetime, fn) && ~isempty(m.lifetime.(fn))
            r.lifetime.(fn) = [r.lifetime.(fn); m.lifetime.(fn)(:)];
            r.lifetime.All = [r.lifetime.All; m.lifetime.(fn)(:)];
        end
    end
    
    % Track group counts
    if isfield(m, 'groupCounts') && ~isempty(m.groupCounts)
        r.groupCounts = [r.groupCounts; m.groupCounts(:)'];
    end
end

function printDataSummary(results)
    fprintf('\n--- Data Summary ---\n');
    fprintf('Control:\n');
    fprintf('  Assembly G1: n=%d\n', length(results.control.assembly.G1));
    fprintf('  Assembly G4: n=%d\n', length(results.control.assembly.G4));
    fprintf('  Assembly G8: n=%d\n', length(results.control.assembly.G8));
    fprintf('  Lifetime G1: n=%d\n', length(results.control.lifetime.G1));
    fprintf('  Lifetime G2: n=%d\n', length(results.control.lifetime.G2));
    fprintf('  Lifetime FA: n=%d\n', length(results.control.lifetime.totalFA));
    
    fprintf('MEKi:\n');
    fprintf('  Assembly G1: n=%d\n', length(results.meki.assembly.G1));
    fprintf('  Assembly G4: n=%d\n', length(results.meki.assembly.G4));
    fprintf('  Assembly G8: n=%d\n', length(results.meki.assembly.G8));
    fprintf('  Lifetime G1: n=%d\n', length(results.meki.lifetime.G1));
    fprintf('  Lifetime G2: n=%d\n', length(results.meki.lifetime.G2));
    fprintf('  Lifetime FA: n=%d\n', length(results.meki.lifetime.totalFA));
end

function generatePlots(results, outputDir)
    fig = figure('Position', [100 100 1500 500], 'Color', 'white');
    cControl = [0.3 0.4 0.7]; cMEKi = [0.8 0.3 0.3];
    
    % Panel D: Assembly (G1, G4, G8, All)
    subplot(1,3,1);
    plotData({results.control.assembly.G1, results.meki.assembly.G1, ...
              results.control.assembly.G4, results.meki.assembly.G4, ...
              results.control.assembly.G8, results.meki.assembly.G8, ...
              results.control.assembly.All, results.meki.assembly.All}, ...
             {'NA Turnover', 'NAs Stalling', 'FAs Stable', 'All Groups'}, cControl, cMEKi);
    ylabel('Assembly Rate (a.u./min)', 'FontWeight', 'bold');
    title('D. Assembly Rates', 'FontSize', 12, 'FontWeight', 'bold');
    
    % Add sample sizes as text
    addSampleSizes(results.control.assembly, results.meki.assembly, {'G1', 'G4', 'G8', 'All'});
    
    % Panel E: Disassembly (G1, G4, G8, All)
    subplot(1,3,2);
    plotData({results.control.disassembly.G1, results.meki.disassembly.G1, ...
              results.control.disassembly.G4, results.meki.disassembly.G4, ...
              results.control.disassembly.G8, results.meki.disassembly.G8, ...
              results.control.disassembly.All, results.meki.disassembly.All}, ...
             {'NA Turnover', 'NAs Stalling', 'FAs Stable', 'All Groups'}, cControl, cMEKi);
    ylabel('Disassembly Rate (a.u./min)', 'FontWeight', 'bold');
    title('E. Disassembly Rates', 'FontSize', 12, 'FontWeight', 'bold');
    
    % Add sample sizes as text
    addSampleSizes(results.control.disassembly, results.meki.disassembly, {'G1', 'G4', 'G8', 'All'});
    
    % Panel F: Lifetime (G1, G2, totalFA, All)
    subplot(1,3,3);
    plotData({results.control.lifetime.G1, results.meki.lifetime.G1, ...
              results.control.lifetime.G2, results.meki.lifetime.G2, ...
              results.control.lifetime.totalFA, results.meki.lifetime.totalFA, ...
              results.control.lifetime.All, results.meki.lifetime.All}, ...
             {'NA Failing', 'NA Maturing', 'Total FA', 'All Groups'}, cControl, cMEKi);
    ylabel('Lifetime (sec)', 'FontWeight', 'bold');
    title('F. Adhesion Lifetime', 'FontSize', 12, 'FontWeight', 'bold');
    
    % Add sample sizes as text
    addSampleSizes(results.control.lifetime, results.meki.lifetime, {'G1', 'G2', 'totalFA', 'All'});
    
    % Add legend
    legend({'Control', 'MEKi'}, 'Location', 'northeast');
    
    % Add overall title
    sgtitle('Control vs MEKi Focal Adhesion Comparison (* p<0.05, ** p<0.01, *** p<0.001)', ...
            'FontSize', 14, 'FontWeight', 'bold');
    
    % Save in multiple formats
    saveas(fig, fullfile(outputDir, 'Comparison.png'));
    saveas(fig, fullfile(outputDir, 'Comparison.fig'));
    
    % Also save as PDF for publication
    set(fig, 'PaperPositionMode', 'auto');
    print(fig, fullfile(outputDir, 'Comparison.pdf'), '-dpdf', '-r300');
    
    % Save statistical summary
    saveStatsSummary(results, outputDir);
    
    fprintf('Plots saved:\n');
    fprintf('  - %s\n', fullfile(outputDir, 'Comparison.png'));
    fprintf('  - %s\n', fullfile(outputDir, 'Comparison.pdf'));
    fprintf('  - %s\n', fullfile(outputDir, 'StatsSummary.txt'));
end

function addSampleSizes(controlData, mekiData, fields)
    % Add sample size annotations at bottom of plot
    ax = gca;
    yLim = ylim;
    yPos = yLim(1) - (yLim(2) - yLim(1)) * 0.08;
    
    for g = 1:length(fields)
        fn = fields{g};
        nC = length(controlData.(fn));
        nM = length(mekiData.(fn));
        text(g, yPos, sprintf('n=%d/%d', nC, nM), ...
            'HorizontalAlignment', 'center', 'FontSize', 7, 'Color', [0.4 0.4 0.4]);
    end
end

function saveStatsSummary(results, outputDir)
    % Save statistical summary to text file
    fid = fopen(fullfile(outputDir, 'StatsSummary.txt'), 'w');
    
    fprintf(fid, 'FOCAL ADHESION ANALYSIS - STATISTICAL SUMMARY\n');
    fprintf(fid, 'Generated: %s\n\n', datestr(now));
    
    % Assembly
    fprintf(fid, '=== ASSEMBLY RATES ===\n');
    printStats(fid, 'G1 (NA Turnover)', results.control.assembly.G1, results.meki.assembly.G1);
    printStats(fid, 'G4 (Stalling)', results.control.assembly.G4, results.meki.assembly.G4);
    printStats(fid, 'G8 (Stable FA)', results.control.assembly.G8, results.meki.assembly.G8);
    printStats(fid, 'All Groups', results.control.assembly.All, results.meki.assembly.All);
    
    % Disassembly
    fprintf(fid, '\n=== DISASSEMBLY RATES ===\n');
    printStats(fid, 'G1 (NA Turnover)', results.control.disassembly.G1, results.meki.disassembly.G1);
    printStats(fid, 'G4 (Stalling)', results.control.disassembly.G4, results.meki.disassembly.G4);
    printStats(fid, 'G8 (Stable FA)', results.control.disassembly.G8, results.meki.disassembly.G8);
    printStats(fid, 'All Groups', results.control.disassembly.All, results.meki.disassembly.All);
    
    % Lifetime
    fprintf(fid, '\n=== LIFETIME ===\n');
    printStats(fid, 'G1 (NA Failing)', results.control.lifetime.G1, results.meki.lifetime.G1);
    printStats(fid, 'G2 (NA Maturing)', results.control.lifetime.G2, results.meki.lifetime.G2);
    printStats(fid, 'Total FA', results.control.lifetime.totalFA, results.meki.lifetime.totalFA);
    printStats(fid, 'All Groups', results.control.lifetime.All, results.meki.lifetime.All);
    
    fclose(fid);
end

function printStats(fid, name, ctrl, meki)
    fprintf(fid, '\n%s:\n', name);
    if ~isempty(ctrl)
        fprintf(fid, '  Control: %.3f ± %.3f (n=%d)\n', mean(ctrl), std(ctrl), length(ctrl));
    else
        fprintf(fid, '  Control: No data\n');
    end
    if ~isempty(meki)
        fprintf(fid, '  MEKi: %.3f ± %.3f (n=%d)\n', mean(meki), std(meki), length(meki));
    else
        fprintf(fid, '  MEKi: No data\n');
    end
    
    % Statistical test
    if ~isempty(ctrl) && ~isempty(meki) && length(ctrl) > 1 && length(meki) > 1
        [~, p] = ttest2(ctrl, meki);
        fprintf(fid, '  t-test p-value: %.4f', p);
        if p < 0.001
            fprintf(fid, ' (***)\n');
        elseif p < 0.01
            fprintf(fid, ' (**)\n');
        elseif p < 0.05
            fprintf(fid, ' (*)\n');
        else
            fprintf(fid, ' (ns)\n');
        end
    end
end

function plotData(data, labels, c1, c2)
    hold on;
    yMaxAll = 0;  % Track overall maximum for significance brackets
    
    for g = 1:length(labels)
        d1 = data{(g-1)*2+1}; d2 = data{(g-1)*2+2};
        
        % Remove NaN and Inf values
        d1 = d1(isfinite(d1));
        d2 = d2(isfinite(d2));
        
        p1 = g-0.2; p2 = g+0.2;
        
        if ~isempty(d1)
            scatter(p1 + (rand(length(d1),1)-0.5)*0.15, d1, 15, c1, 'filled', 'MarkerFaceAlpha', 0.5);
            if length(d1) > 1, drawBox(d1, p1, c1); end
            yMaxAll = max(yMaxAll, max(d1));
        end
        if ~isempty(d2)
            scatter(p2 + (rand(length(d2),1)-0.5)*0.15, d2, 15, c2, 'filled', 'MarkerFaceAlpha', 0.5);
            if length(d2) > 1, drawBox(d2, p2, c2); end
            yMaxAll = max(yMaxAll, max(d2));
        end
    end
    
    % Add significance brackets AFTER all data is plotted
    for g = 1:length(labels)
        d1 = data{(g-1)*2+1}; d2 = data{(g-1)*2+2};
        
        % Remove NaN and Inf values
        d1 = d1(isfinite(d1));
        d2 = d2(isfinite(d2));
        
        p1 = g-0.2; p2 = g+0.2;
        
        if ~isempty(d1) && ~isempty(d2) && length(d1) > 1 && length(d2) > 1
            % Perform t-test
            [~, p] = ttest2(d1, d2);
            
            % Determine significance level
            if p < 0.001
                sig = '***';
            elseif p < 0.01
                sig = '**';
            elseif p < 0.05
                sig = '*';
            else
                sig = 'ns';
            end
            
            % Draw significance bracket
            localMax = max([max(d1), max(d2)]);
            bracketY = localMax * 1.1;
            bracketTop = localMax * 1.15;
            
            % Draw bracket lines
            line([p1, p1], [bracketY, bracketTop], 'Color', 'k', 'LineWidth', 1);
            line([p2, p2], [bracketY, bracketTop], 'Color', 'k', 'LineWidth', 1);
            line([p1, p2], [bracketTop, bracketTop], 'Color', 'k', 'LineWidth', 1);
            
            % Add significance text
            if strcmp(sig, 'ns')
                text(g, bracketTop * 1.02, sig, 'HorizontalAlignment', 'center', ...
                    'FontSize', 9, 'Color', [0.5 0.5 0.5]);
            else
                text(g, bracketTop * 1.02, sig, 'HorizontalAlignment', 'center', ...
                    'FontSize', 14, 'FontWeight', 'bold');
            end
            
            % Add p-value below (smaller font)
            if p < 0.001
                pStr = 'p<0.001';
            else
                pStr = sprintf('p=%.3f', p);
            end
            text(g, bracketTop * 1.08, pStr, 'HorizontalAlignment', 'center', ...
                'FontSize', 7, 'Color', [0.4 0.4 0.4]);
        end
    end
    
    set(gca, 'XTick', 1:length(labels), 'XTickLabel', labels);
    xlim([0.5 length(labels)+0.5]); 
    
    % Adjust y-axis to fit brackets - handle case where yMaxAll is still 0
    if yMaxAll > 0
        ylim([0, yMaxAll * 1.25]);
    else
        ylim([0, 1]); % Default if no data
    end
    
    box off; 
    grid on;
end

function drawBox(d, pos, c)
    q = quantile(d, [0.25 0.5 0.75]);
    iqr = q(3) - q(1);
    
    % Box (IQR)
    rectangle('Position', [pos-0.1, q(1), 0.2, iqr], 'EdgeColor', c, 'LineWidth', 1.5);
    
    % Median line
    line([pos-0.1, pos+0.1], [q(2), q(2)], 'Color', c, 'LineWidth', 2.5);
    
    % Whiskers (to min/max within 1.5*IQR)
    lowerWhisker = max(min(d), q(1) - 1.5*iqr);
    upperWhisker = min(max(d), q(3) + 1.5*iqr);
    
    % Lower whisker
    line([pos, pos], [q(1), lowerWhisker], 'Color', c, 'LineWidth', 1);
    line([pos-0.05, pos+0.05], [lowerWhisker, lowerWhisker], 'Color', c, 'LineWidth', 1);
    
    % Upper whisker  
    line([pos, pos], [q(3), upperWhisker], 'Color', c, 'LineWidth', 1);
    line([pos-0.05, pos+0.05], [upperWhisker, upperWhisker], 'Color', c, 'LineWidth', 1);
end