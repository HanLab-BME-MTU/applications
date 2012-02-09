function exitflag = algorithmIncremental(dirPath,rootPathWin,rootPathUnix)
% Disable stupid MATLAB parfor warning!
warning('off','MATLAB:mir_warning_maybe_uninitialized_temporary');

tStart = tic;

smallClusterSize = 3;

% INIT DATA
disp('-- INIT ----------------------------------------');

dirPath = checkPath(dirPath,rootPathWin,rootPathUnix);
list = dir(dirPath);

% Find the first configuration file in the directory
for i=3:numel(list)
    if strcmp(list(i).name(end-3:end),'.cfg')
        configName = list(i).name;
        break;
    else
        exitflag = 'Main: Configuration file not found!';
        return;
    end
end

% Load the configuration file
cfg = Config.load([dirPath configName]);
cfg.path = checkPath(cfg.path,rootPathWin,rootPathUnix);

if isunix
    disp('Main: Unix system detected: Display and snapshots are disabled!');
    cfg.displayEnabled = false;
    cfg.snapshotsEnabled = false;
end

% Initialize Imaris
im = Imaris(cfg.displayEnabled,cfg.snapshotsEnabled,cfg.snapshotsPath);
im.setupScene();

% Check if *.i.dat file exists
dataName = [];
for i=3:numel(list)
    if strcmp(list(i).name(end-5:end),'.i.dat')
        dataName = list(i).name;
        break;
    end
end

if isempty(dataName) % No *.i.dat file present
    
    % Check if *.d.data file exists
    for i=3:numel(list)
        if strcmp(list(i).name(end-5:end),'.d.dat')
            dataName = list(i).name;
            break;
        end
    end
    
    if isempty(dataName) % No *.d.dat file present
        % Read data
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
        
        % Prefilter data
        if cfg.dataReductionEnabled
            pro.dataReduction(cfg.reductionEdgeRadius,cfg.nReductionRun);
            dis.points();
            im.takeSnapshotAndResetScene(); % ========== SNAPSHOT ============
        end
        
        if cfg.densityFilteringEnabled
            pro.densityFilter(cfg.nNeighborsThreshold,cfg.neighborBallRadius);
            dis.points();
            im.takeSnapshotAndResetScene(); % ========== SNAPSHOT ============
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
    pro.computeOrientation(cfg.filterLength,cfg.angularSampling);
    dis.orientation(100);
    
    im.takeSnapshotAndResetScene(); % ========== SNAPSHOT ============
    
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

pro.initEdges(cfg.initialEdgeRadius);
pro.dissolveClustersSmallerThan(smallClusterSize);
pro.initModels();

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

    pro.assignPointsToModels(cfg.nSigmaThreshold,cfg.modeVar);
    pro.dissolveClustersSmallerThan(smallClusterSize);
    pro.updateModels(cfg.maxCurvature,cfg.fitMethod,cfg.modeVar);

    fprintf('Main: Relaxed %d points in run %d: \n',nnz(parentsOld ~= data.parents),relax);
    
    if cfg.snapshotsEnabled
        dis.points(); dis.models(); % dis.clusters();
    end
    
    im.takeSnapshotAndResetScene(); % ========== SNAPSHOT ============
end

% for curveDegree=1:cfg.maxDegreeBezier+1
for curveDegree=1:cfg.maxDegreeBezier
    fprintf('Main: Current maximal curve degree: %u\n',curveDegree);
    for mergeIter=1:cfg.maxIterMerge
        pro.updateEdges();
        fprintf('Main: Merge run: %d\n',mergeIter);
        
        if curveDegree == cfg.maxDegreeBezier+1
            if pro.mergeClustersBIC(cfg.maxDegreeBezier,-cfg.maxCurvature,cfg.fitMethod,cfg.betaVar,cfg.modeVar)
                disp('Main: Nothing to merge!');
                break;
            end
        else
            if pro.mergeClustersBIC(curveDegree,cfg.maxCurvature,cfg.fitMethod,cfg.betaVar,cfg.modeVar)
            disp('Main: Nothing to merge!');
            break;
            end
        end
        
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
            
            pro.assignPointsToModels(cfg.nSigmaThreshold,cfg.modeVar);
            pro.dissolveClustersSmallerThan(smallClusterSize);
            pro.updateModels(cfg.maxCurvature,cfg.fitMethod,cfg.modeVar);
            
            fprintf('Main: Relaxed %d points in run %d: \n',nnz(parentsOld ~= data.parents),relax);
            
            if cfg.snapshotsEnabled
                dis.points(); dis.models(); % dis.clusters();
            end
            
            im.takeSnapshotAndResetScene(); % ========== SNAPSHOT ============
        end
        
        if mergeIter == cfg.maxIterMerge
            disp('Main: maxIterMerge reached!')
            tEnd = round(toc(tStart));
            data.runTime = tEnd;
            data.save([dirPath cfg.configName '.p.dat']);
            exitflag = 'maxIterMerge reached!';
            return;
        end
        
    end
end

tEnd = round(toc(tStart));
data.runTime = tEnd;
data.save([dirPath cfg.configName '.p.dat']);
exitflag = sprintf('Success! %s',secs2hms(tEnd));

disp('------------------------------------------------');

end




