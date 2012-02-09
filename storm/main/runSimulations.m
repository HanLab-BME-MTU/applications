
%% GENERATE CONFIG

clear all; clc;

nCurves = 100; % 100
nPoints = 40;
offsetY = 500;

types = 4:6;
for type=types
    switch(type)
        case 1
            VARS = 1; % Line
        case 2
            VARS = 1; % Arc
        case 3
            VARS = 1; % Spline
        case 4
            VARS = [10,20,50,100,200]; % Noise
        case 5
            VARS = [5,10,15,20,30,40,50]; % Distance
        case 6
            VARS = [15,30,45,60,75,90]; % Angle
    end
        
    % Get a list with all the .cfg-files in _sim
    path = 'Y:\fsm\harvard\data\Zhuang\_sim\';
    list = dir(path);
    itemIdx = 0;
    nConfig = numel(list)-2;
    items = cell(nConfig,2);
    for i=3:numel(list)
        if ~list(i).isdir
            if strcmp(list(i).name(end-3:end),'.cfg')
                % Config found
                itemIdx = itemIdx + 1;
                items{itemIdx,1} = [path list(i).name];
                items{itemIdx,2} = list(i).name(1:end-4);
            end
        end
    end
    items = items(1:itemIdx,:);
    
    for VAR=VARS
        
        for k=1:size(items,1)
            % Read the main cfg file
            configPath = items{k,1};
            cfg = Config.load(configPath);

            % Generate the data
            dat = Data();
            dis = Display(dat);
            sim = Simulation(dat);

            % sim.setSamplingToRandom();
            sim.setSamplingToRegular();
            
            for c=1:nCurves
                switch(type)
                    case 1 % Line
                        cP = [0 0 0; 200 0 0];
                    case 2 % Arc
                        % cP = [0 0 0; 100 50 0; 200 0 0];
                        % cP = [0 0 0; 100 100 0; 200 0 0];
                        cP = [0 0 0; 100 75 0; 200 0 0];
                    case 3 % Spline
                        cP = [0 0 0; 100 100 0; 200 -100 0; 300 0 0];
                        % cP = [0 0 0; 100 150 0; 200 -150 0; 300 0 0];
                    case 4 % Noise
                        cP = [0 0 0; 100 100 0; 200 -100 0; 300 0 0];
                        sim.setDomain([-100 -200+(c-1)*offsetY 0],[500 400 0]);
                        sim.addRandomNoise(VAR);
                    case 5 % Distance
                        % cP1 = [0 0 0; 100 100 0; 200 -100 0; 300 0 0]; % Spline
                        cP1 = [0 0 0; 200 0 0]; % Line
                        cP2 = cP1;
                        cP1(:,2) = cP1(:,2)+VAR/2;
                        cP2(:,2) = cP2(:,2)-VAR/2;
                        cP1(:,2) = cP1(:,2)+(c-1)*offsetY;
                        sim.bezier(cP1,nPoints);
                        % sim.bezierWithOrthogonalGaussianNoise2D(cP1,nPoints,cfg.modeVar);
                        dat.simModelBezCP = [dat.simModelBezCP {cP1}];
                        cP = cP2;
                    case 6 % Angle
                        R_plus = [cosd(VAR/2) -sind(VAR/2); sind(VAR/2) cosd(VAR/2)];
                        R_min = [cosd(-VAR/2) -sind(-VAR/2); sind(-VAR/2) cosd(-VAR/2)];
                        
                        cP1 = [0 0 0; 200 0 0];
                        cP2 = cP1;
                        
                        center = (cP1(1,1:2)+cP1(end,1:2))/2;
                        cP1(:,1:2) = cP1(:,1:2)-[center;center];
                        cP2(:,1:2) = cP2(:,1:2)-[center;center];
                        
                        cP1(:,1:2) = cP1(:,1:2)*R_plus;
                        cP2(:,1:2) = cP2(:,1:2)*R_min;
                        
                        cP1(:,1:2) = cP1(:,1:2)+[center;center];
                        cP2(:,1:2) = cP2(:,1:2)+[center;center];
                        
                        cP1(:,2) = cP1(:,2)+(c-1)*offsetY;
                        sim.bezier(cP1,nPoints);
                        % sim.bezierWithOrthogonalGaussianNoise2D(cP1,nPoints,cfg.modeVar);
                        dat.simModelBezCP = [dat.simModelBezCP {cP1}];
                        cP = cP2;
                end
                cP(:,2) = cP(:,2)+(c-1)*offsetY;
                sim.bezier(cP,nPoints);
                % sim.bezierWithOrthogonalGaussianNoise2D(cP,nPoints,cfg.modeVar);
                dat.simModelBezCP = [dat.simModelBezCP {cP}];
            end
            
            sim.addGaussianNoise([cfg.errorX cfg.errorY 0]);
            pro = Processor(dat);
            pro.setErrorArray(cfg.errorX,cfg.errorY,cfg.errorZ);
            
            % Display data
            dis.points();
            
            % Create folder and save cfg and dat file
            switch(type)
                case 1
                    configName = ['line_' items{k,2}];
                case 2
                    configName = ['arc_' items{k,2}];
                case 3
                    configName = ['spline_' items{k,2}];
                case 4
                    configName = [sprintf('noise_%04u_',VAR) items{k,2}];
                case 5
                    configName = [sprintf('distance_%04u_',VAR) items{k,2}];
                case 6
                    configName = [sprintf('angle_%04u_',VAR) items{k,2}];
            end
            cfg.configName = configName;
            parentFolder = 'Y:\fsm\harvard\data\Zhuang\_sim\';
            folderName = [configName '-'];
            mkdir(parentFolder,folderName);
            configPathOut = ['Y:\fsm\harvard\data\Zhuang\_sim\' folderName '\' configName '.cfg'];
            cfg.save(configPathOut);
            dataPathOut = ['Y:\fsm\harvard\data\Zhuang\_sim\' folderName '\' configName '.d.dat'];
            dat.save(dataPathOut);
            
        end
    end
    
    disp('==========================')
end


%% EVALUATE DATA

clear all; clc; format short;

offsetY = 500;

% Get a list with all the .dat-files in _sim
path = 'Y:\fsm\harvard\data\Zhuang\_sim\';
list = dir(path);
itemIdx = 0;
nDataSets = numel(list)-2;
items = cell(nDataSets,2);
for i=3:numel(list)
    if list(i).isdir
        sublist = dir([path list(i).name '\']);
        for k=3:numel(sublist)
            if strcmp(sublist(k).name(end-3:end),'.cfg')
                % Config found - look for corresponding data file
                for j=3:numel(sublist)
                    if strcmp(sublist(j).name(end-5:end),'.p.dat')
                        itemIdx = itemIdx + 1;
                        items{itemIdx,1} = [path list(i).name '\' sublist(k).name];
                        items{itemIdx,2} = [path list(i).name '\' sublist(j).name];
                        items{itemIdx,3} = sublist(j).name;
                        break;
                    end
                end
                break;
            end
        end
    end
end
items = items(1:itemIdx,:);

output = cell(4,size(items,1)); 

% Loop through all the files
lineCounter = 1;    
arcCounter = 1;   
splineCounter = 1;  
distanceCounter = 1;  
angleCounter = 1;  
noiseCounter = 1;  
for i=1:size(items,1)
    
    % Read the .dat-file and .cfg-file
    cfg = Config.load(items{i,1});
    dat = Data.load(items{i,2});
    pro = Processor(dat);
    pro.dissolveClustersSmallerThan(5);
    
    fprintf('-------------\nFile: %s\n',items{i,1});
    
    modelBezCP = dat.modelBezCP;
    
    nModel = size(modelBezCP,1);
    idxModel = zeros(nModel,1);
    for m=1:nModel
        idxModel(m) = round(mean(modelBezCP{m}(:,2))/offsetY)+1;
    end
    
    nSimModel = length(dat.simModelBezCP);
    idxSimModel = zeros(nSimModel,1);
    for m=1:nSimModel
        idxSimModel(m) = round(mean(dat.simModelBezCP{m}(:,2))/offsetY)+1;
    end
    
    failed = false(max(idxSimModel),1);
    hausdorffDist = zeros(max(idxSimModel),1);
    
    for s=1:max(idxSimModel)
        nSimMod = nnz(idxSimModel == s);
        nMod = nnz(idxModel == s);
        
        if nMod ~= nSimMod
            % The number of models is not the same
            failed(s) = true;
        else
            % Compare the models pairwise
            m=find(idxModel == s);
            n=find(idxSimModel == s);
            
            for k=1:numel(m)
                                
                if size(dat.simModelBezCP{n(k)},1) ~= size(modelBezCP{m(k)},1)
                    % The complexity of the models is not the same
                    % if isempty(strfind(items{i,1},'angle')) && isempty(strfind(items{i,1},'distance')) % All except angle and distance
                    if isempty(strfind(items{i,1},'angle')) % All except angle
                        failed(s) = true;
                    end
                end
            end
            if ~isempty(strfind(items{i,1},'angle')) % Angle only
                cP1 = modelBezCP{m(1)};
                cP2 = modelBezCP{m(2)};
                dist = segments_dist_3d (cP1(1,:)',cP1(end,:)',cP2(1,:)',cP2(end,:)');
                if dist > 1
                    failed(s) = true;
                end
            end
        end
                
        if failed(s) == false
            % Compute the model distance
            nSamples = 1000;
            refPointsCell = cellfun(@(a) renderBezier(a,linspace(0,1,nSamples)'),dat.simModelBezCP(idxSimModel == s),'UniformOutput',0);
            refPoints = zeros(numel(refPointsCell)*nSamples,3);
            for m=1:numel(refPointsCell)
                refPoints((m-1)*nSamples+1:m*nSamples,:) = refPointsCell{m};
            end
            
            pointsCell = cellfun(@(a) renderBezier(a,linspace(0,1,nSamples)'),modelBezCP(idxModel == s),'UniformOutput',0);
            points = zeros(numel(pointsCell)*nSamples,3);
            for m=1:numel(pointsCell)
                points((m-1)*nSamples+1:m*nSamples,:) = pointsCell{m};
            end
            d = createDistanceMatrix(refPoints,points);
            hausdorffDist(s) = max(min(d));
        end
    end
        
    hausdorffDist = hausdorffDist(~failed);
    
    meanHausdorffDist = mean(hausdorffDist);
    stdDevHausdorffDist = sqrt(1/numel(hausdorffDist)*sum((hausdorffDist-meanHausdorffDist).^2));
    fractionOfFailedModels = nnz(failed)/numel(failed);
    nFailedModels = nnz(failed);
    [meanHausdorffDist;stdDevHausdorffDist;1-fractionOfFailedModels]
    
    %     figure(i);
    %     hist(hausdorffDist,20)
    
    if ~isempty(strfind(items{i,1},'line')) && isempty(strfind(items{i,1},'spline'))
        lineCounter = lineCounter + 1;  
        x = 2; y = lineCounter;
    elseif strfind(items{i,1},'arc')
        arcCounter = arcCounter + 1;  
        x = 8; y = arcCounter;
    elseif strfind(items{i,1},'spline')
        splineCounter = splineCounter + 1;  
        x = 14; y = splineCounter;
    elseif strfind(items{i,1},'distance')
        distanceCounter = distanceCounter + 1;  
        x = 20; y = distanceCounter;
    elseif strfind(items{i,1},'angle')
        angleCounter = angleCounter + 1;
        x = 26; y = angleCounter;
    elseif strfind(items{i,1},'noise')
        noiseCounter = noiseCounter + 1;  
        x = 32; y = noiseCounter;
    end
    % output(x,y) = items(i,3);
    output(x+1,y) = {meanHausdorffDist};
    output(x+2,y) = {stdDevHausdorffDist};
    output(x+3,y) = {1-fractionOfFailedModels};
    
end

xlswrite('C:\Users\PB93\Desktop\output.xls',output);

disp('==========================')








