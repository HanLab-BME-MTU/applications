function exitflag = algorithm(dirPath)
% Disable stupid MATLAB parfor warning!
warning('off','MATLAB:mir_warning_maybe_uninitialized_temporary');

global timerX99; % :-D
timerX99 = Timing();
timerX99.start('Main');

smallClusterSize = 3;
dRef = 40; alpha = 0.25; samplePeriod = 10; dMaxAlong = 160; dMinAway = 20;
 
% INIT DATA
disp('-- INIT ----------------------------------------');
timerX99.start('Read Configuration');

list = dir(dirPath);
list = list(~vertcat(list.isdir));

% Find the first configuration file in the directory
for i=1:numel(list)
    if strcmp(list(i).name(end-3:end),'.cfg')
        configName = list(i).name;
        break;
    end
    if i == numel(list)
        exitflag = 'Main: Configuration file not found!';
        return;
    end
end

% Load the configuration file
cfg = Config.load([dirPath configName]);
cfg.path
cfg.path = checkStormPath(cfg.path);
cfg.path
if isunix
    disp('Main: Unix system detected: Display and snapshots are disabled!');
    cfg.displayEnabled = false;
    cfg.snapshotsEnabled = false;
end
timerX99.stop('Read Configuration');

% Initialize Imaris
im = Imaris(cfg.displayEnabled,cfg.snapshotsEnabled,cfg.snapshotsPath);
im.setupScene();

% Check if *.i.dat file exists
dataName = [];
for i=1:numel(list)
    if strcmp(list(i).name(end-5:end),'.i.dat')
        dataName = list(i).name;
        break;
    end
end

if isempty(dataName) % No *.i.dat file present
    
    % Check if *.d.data file exists
    for i=1:numel(list)
        if strcmp(list(i).name(end-5:end),'.d.dat')
            dataName = list(i).name;
            break;
        end
    end
    
    if isempty(dataName) % No *.d.dat file present
        % Read data
        timerX99.start('Read Data');
        [cfg.path cfg.fileName]
        data = Data.read([cfg.path cfg.fileName]);
        dis = Show(data,im);
        pro = Processor(data);
        pro.cropData(cfg.roiPosition,cfg.roiSize);
        
        % Check if region is empty
        if data.nPoints < 10
            exitflag = 'Main: Not enough data points!';
            return;
        end
        
        pro.centerData();
        dis.points();
        im.fitAndSaveCamera();
        im.takeSnapshotAndResetScene(); % ========== SNAPSHOT ============
        data.rawPoints = data.points;
        timerX99.stop('Read Data');
        
        % Prefilter data
        if cfg.dataReductionEnabled
            timerX99.start('Data Reduction');
            pro.dataReduction(cfg.reductionEdgeRadius,cfg.nReductionRun);
            dis.points();
            im.takeSnapshotAndResetScene(); % ========== SNAPSHOT ============
            timerX99.stop('Data Reduction');
        end
        
        if cfg.densityFilteringEnabled
            timerX99.start('Density Filter');
            pro.densityFilter(cfg.nNeighborsThreshold,cfg.neighborBallRadius);
            dis.points();
            im.takeSnapshotAndResetScene(); % ========== SNAPSHOT ============
            timerX99.stop('Density Filter');
        end
        
        pro.setErrorArray(cfg.errorX,cfg.errorY,cfg.errorZ);
        data.save([dirPath cfg.configName '.d.dat']);
        
    else % Load *.d.dat file
        data = Data.load([dirPath dataName]);
        dis = Show(data,im);
        pro = Processor(data);
        dis.points();
        im.fitAndSaveCamera();
        im.takeSnapshotAndResetScene(); % ========== SNAPSHOT ============
    end
    
    disp('------------------------------------------------');
    disp('-- PRE-PROCESS ---------------------------------');
        
    dis.points();
    timerX99.start('Orientation');
    pro.computeOrientation(cfg.filterLength,cfg.angularSampling);
    timerX99.stop('Orientation');
    dis.orientation(100);
        
    im.takeSnapshotAndResetScene(); % ========== SNAPSHOT ============
    
    timerX99.start('Geometric Matching');
    pro.initClusters();
    pro.initEdges(cfg.initialEdgeRadiusGeom);
    dis.initialEdges();
    
    for k=1:cfg.nIterGeomMatching
        pro.mergeClustersGeom(cfg.modelLength,cfg.angleThreshold);
        pro.saveEdgesToHistory();
        pro.saveClustersToHistory();
        % dis.edgesHistoryLast();
        pro.updateEdges();
        dis.points();
        im.takeSnapshotAndResetScene(); % ========== SNAPSHOT ============
    end
    timerX99.stop('Geometric Matching');

    data.save([dirPath cfg.configName '.i.dat']);
    
else % Load *.i.dat file
    
    data = Data.load([dirPath dataName]);
    dis = Show(data,im);
    pro = Processor(data);
    dis.points();
    im.fitAndSaveCamera();
    im.takeSnapshotAndResetScene(); % ========== SNAPSHOT ============
    
end

if cfg.edgeWidthInitFree > 0
    pro.dissolveClustersCloseToEdge(cfg.edgeWidthInitFree);
    % idxPointsEdge = pro.dissolveClustersCloseToEdge(cfg.edgeWidthInitFree);
    % dis.imaris.displayPoints(data.points(idxPointsEdge,:),8,[1.0 0.0 0.0 0.0],'Data Set Edge');
end

disp('------------------------------------------------');
disp('-- PROCESS -------------------------------------');

timerX99.start('Process');
timerX99.start('Init Clusters');
pro.initEdges(cfg.initialEdgeRadius);
pro.dissolveClustersSmallerThan(smallClusterSize);

timerX99.start('Init Models');
pro.initModels(cfg.betaVar,cfg.modeVar);
timerX99.stop('Init Models');

if cfg.snapshotsEnabled
    dis.points(); dis.models(); % dis.clusters();
end

im.takeSnapshotAndResetScene(); % ========== SNAPSHOT ============

parentsOld = zeros(size(data.parents));
relax = 0;
while any(parentsOld ~= data.parents)
    relax = relax+1;
    if relax > cfg.nIterEM
        break;
    end
    parentsOld = data.parents;

    pro.assignPointsToModels(cfg.nSigmaThreshold);

    pro.dissolveClustersSmallerThan(smallClusterSize);
    pro.updateModels(cfg.maxCurvature,cfg.fitMethod,cfg.betaVar,cfg.modeVar);

    fprintf('Main: Relaxed %d points in run %d: \n',nnz(parentsOld ~= data.parents),relax);

    if cfg.snapshotsEnabled
        dis.points(); dis.models(); % dis.clusters();
    end

    im.takeSnapshotAndResetScene(); % ========== SNAPSHOT ============
end
timerX99.stop('Init Clusters');
timerX99.start('Main Loop');
for mergeIter=1:cfg.maxIterMerge
    
%     tic
%     pro.updateEdges();
%     b = size(data.edges,1)
%     toc
    
    tic
    pro.updateEdgesAnisotropic(dRef,alpha,samplePeriod,dMaxAlong,dMinAway);
    c = size(data.edges,1);
    toc
    
%     tic
%     pro.updateEdgesEndPoints(cfg.initialEdgeRadius);
%     b = size(data.edges,1)
%     toc
    
    fprintf('Main: Merge run: %d\n',mergeIter);
    
    timerX99.start('Compute BIC Weights');
    if pro.mergeClustersBIC(cfg.maxDegreeBezier,cfg.maxCurvature,cfg.fitMethod,cfg.betaVar,cfg.modeVar)
        disp('Main: Nothing to merge!');
        timerX99.stop('Compute BIC Weights',data.nClusters);
        break;
    else
        timerX99.stop('Compute BIC Weights',data.nClusters);
    end

    if cfg.snapshotsEnabled
        dis.points(); dis.models(); % dis.clusters();
    end

    im.takeSnapshotAndResetScene(); % ========== SNAPSHOT ============

    parentsOld = zeros(size(data.parents));
    relax = 0;
    timerX99.start('Refine Model');
    while any(parentsOld ~= data.parents)
        relax = relax+1;
        if relax > cfg.nIterEM
            break;
        end
        parentsOld = data.parents;

        timerX99.start('Assign Points');
        pro.assignPointsToModels(cfg.nSigmaThreshold);
        pro.dissolveClustersSmallerThan(smallClusterSize);
        timerX99.stop('Assign Points');
    
        timerX99.start('Update Models');
        pro.updateModels(cfg.maxCurvature,cfg.fitMethod,cfg.betaVar,cfg.modeVar);
        timerX99.stop('Update Models');

        fprintf('Main: Relaxed %d points in run %d: \n',nnz(parentsOld ~= data.parents),relax);

        if cfg.snapshotsEnabled
            dis.points(); dis.models(); % dis.clusters();
        end

        im.takeSnapshotAndResetScene(); % ========== SNAPSHOT ============
    end
    timerX99.stop('Refine Model');

    if mergeIter == cfg.maxIterMerge
        disp('Main: maxIterMerge reached!')
        data.runTime = timerX99.stop('Main');
        data.save([dirPath cfg.configName '.p.dat']);
        exitflag = 'maxIterMerge reached!';
        return;
    end
    
    timerX99.save([dirPath cfg.configName '_tmp.tim']);

end
timerX99.stop('Main Loop');
timerX99.stop('Process');

data.runTime = timerX99.stop('Main');
timerX99.running()
timerX99.save([dirPath cfg.configName '.tim']);
data.save([dirPath cfg.configName '.p.dat']);
exitflag = ['Success! ' cell2mat(timerX99.getFormatted('Main'))];

disp('------------------------------------------------');

end




