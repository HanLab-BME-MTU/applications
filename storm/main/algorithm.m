function exitflag = algorithm(dirPath)
% Disable stupid MATLAB parfor warning!
warning('off','MATLAB:mir_warning_maybe_uninitialized_temporary');

global stormTimer__; % :-D
stormTimer__ = Timing();
stormTimer__.start('Main');

smallClusterSize = 3;
dRef = 40; alpha = 0.25; samplePeriod = 10; dMaxAlong = 160; dMinAway = 20;
 
% INIT DATA
disp('-- INIT ----------------------------------------');
stormTimer__.start('Read Configuration');

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
pos = strfind(cfg.path,'_data');
cfg.path = [getStormPath() strrep(strrep(cfg.path(pos:end),'\',filesep),'/',filesep)];

if isunix
    disp('Main: Unix system detected: Display and snapshots are disabled!');
    cfg.displayEnabled = false;
    cfg.snapshotsEnabled = false;
end
stormTimer__.stop('Read Configuration');

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
        stormTimer__.start('Read Data');
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
        stormTimer__.stop('Read Data');
        
        % Prefilter data
        if cfg.dataReductionEnabled
            stormTimer__.start('Data Reduction');
            pro.dataReduction(cfg.reductionEdgeRadius,cfg.nReductionRun);
            dis.points();
            im.takeSnapshotAndResetScene(); % ========== SNAPSHOT ============
            stormTimer__.stop('Data Reduction');
        end
        
        if cfg.densityFilteringEnabled
            stormTimer__.start('Density Filter');
            pro.densityFilter(cfg.nNeighborsThreshold,cfg.neighborBallRadius);
            dis.points();
            im.takeSnapshotAndResetScene(); % ========== SNAPSHOT ============
            stormTimer__.stop('Density Filter');
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
    stormTimer__.start('Orientation');
    pro.computeOrientation(cfg.filterLength,cfg.angularSampling);
    stormTimer__.stop('Orientation');
    dis.orientation(100);
        
    im.takeSnapshotAndResetScene(); % ========== SNAPSHOT ============
    
    stormTimer__.start('Geometric Matching');
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
    stormTimer__.stop('Geometric Matching');

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

stormTimer__.start('Process');
stormTimer__.start('Init Clusters');
pro.initEdges(cfg.initialEdgeRadius);
pro.dissolveClustersSmallerThan(smallClusterSize);

stormTimer__.start('Init Models');
pro.initModels(cfg.betaVar,cfg.modeVar);
stormTimer__.stop('Init Models');

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
stormTimer__.stop('Init Clusters');
stormTimer__.start('Main Loop');
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
    
    stormTimer__.start('Compute BIC Weights');
    if pro.mergeClustersBIC(cfg.maxDegreeBezier,cfg.maxCurvature,cfg.fitMethod,cfg.betaVar,cfg.modeVar)
        disp('Main: Nothing to merge!');
        stormTimer__.stop('Compute BIC Weights',data.nClusters);
        break;
    else
        stormTimer__.stop('Compute BIC Weights',data.nClusters);
    end

    if cfg.snapshotsEnabled
        dis.points(); dis.models(); % dis.clusters();
    end

    im.takeSnapshotAndResetScene(); % ========== SNAPSHOT ============

    parentsOld = zeros(size(data.parents));
    relax = 0;
    stormTimer__.start('Refine Model');
    while any(parentsOld ~= data.parents)
        relax = relax+1;
        if relax > cfg.nIterEM
            break;
        end
        parentsOld = data.parents;

        stormTimer__.start('Assign Points');
        pro.assignPointsToModels(cfg.nSigmaThreshold);
        pro.dissolveClustersSmallerThan(smallClusterSize);
        stormTimer__.stop('Assign Points');
    
        stormTimer__.start('Update Models');
        pro.updateModels(cfg.maxCurvature,cfg.fitMethod,cfg.betaVar,cfg.modeVar);
        stormTimer__.stop('Update Models');

        fprintf('Main: Relaxed %d points in run %d: \n',nnz(parentsOld ~= data.parents),relax);

        if cfg.snapshotsEnabled
            dis.points(); dis.models(); % dis.clusters();
        end

        im.takeSnapshotAndResetScene(); % ========== SNAPSHOT ============
    end
    stormTimer__.stop('Refine Model');

    if mergeIter == cfg.maxIterMerge
        disp('Main: maxIterMerge reached!')
        data.runTime = stormTimer__.stop('Main');
        data.save([dirPath cfg.configName '.p.dat']);
        exitflag = 'maxIterMerge reached!';
        return;
    end
    
    stormTimer__.save([dirPath cfg.configName '_tmp.tim']);
    data.save([dirPath cfg.configName '._tmp.dat']);

end
stormTimer__.stop('Main Loop');
stormTimer__.stop('Process');

data.runTime = stormTimer__.stop('Main');
stormTimer__.running()
stormTimer__.save([dirPath cfg.configName '.tim']);
data.save([dirPath cfg.configName '.p.dat']);
exitflag = ['Success! ' cell2mat(stormTimer__.getFormatted('Main'))];

disp('------------------------------------------------');

end




