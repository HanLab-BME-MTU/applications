%%% Description
% This script shows how to pull tracking results and provides some basic
% plotting of the track output. 

%% %% USER INPUT 
xpName='Metaphase'; % For convenience, Figure window will bear this name.

%% LOADING MOVIES INFO
% MovieList paths: 
MLPath='/project/cellbiology/gdanuser/shared/proudot/project/lattice-track/data-analysis/chemical-inhibition/2015_06_17(EB1_GFP_drugs)/analysis/'

% MovieList FileName (the combination of condition you want to compare). 
movieListFileNames={'controlMetaphase.mat','10nMMetaphase.mat'};
%movieListFileNames={'controlMetaphase.mat','10nMMetaphase.mat','33nMMetaphase.mat','100nMMetaphase.mat'};
%movieListFileNames={'controlInterphase.mat','10nMInterphase.mat','33nMInterphase.mat','100nMInterphase.mat',};
%movieListFileNames={'controlAnaphase.mat','10nMAnaphase.mat','33nMAnaphase.mat'};

% Name associated to each MovieList
conditionName={'control','10nM'};
%conditionName={'control','10nM','33nM','100nM'};

% Build the array of MovieList (auto)
aMovieListArray=[];
CellIdx=cell(1,length(movieListFileNames));
for i=1:length(movieListFileNames)
    aMovieListArray=[aMovieListArray MovieList.loadMatFile([MLPath movieListFileNames{i}])];
    CellIdx{i}=[1:length(aMovieListArray(i).movieDataFile_)];
end;


%% Input parameters
radius=4000;
detectionMethod='pointSourceAutoSigmaLM';
poleScale=3;
poleDetectionMethod=['simplex_scale_' num2str(poleScale,'%03d')];

%% %% MEASUREMENT VECTOR ALLOCATION

% Identifying the maximum number of cell for a single condition accross all condition to help
% results storage
maxCellNb=max(arrayfun(@(x) length(x.movieDataFile_),aMovieListArray)); 

% Per-cell lifetime histograms
lifetimeBins=2:20;
lifetimesHist=cell(length(aMovieListArray),maxCellNb); 

% Per-cell avg track lifetime 
meanLft=nan(length(aMovieListArray),maxCellNb);

% Per-cell avg track length 
meanLength=nan(length(aMovieListArray),maxCellNb);

% Per-cell avg track speed
meanSpeed=nan(length(aMovieListArray),maxCellNb);

% Azimuth and elevation
cumulAzi=cell(length(aMovieListArray),maxCellNb); 
cumulElev=cell(length(aMovieListArray),maxCellNb); 
cumulPoleId=cell(length(aMovieListArray),maxCellNb); 

%conditionName=cell(1,length(aMovieListArray));
%% loop over each different conditions
for k=1:length(aMovieListArray)
    ML=aMovieListArray(k);
    %[~,conditionName{k}]=fileparts(fileparts(ML.movieListPath_));
    %% loop over each different cell in each condition
    parfor i=1:length(ML.movieDataFile_)
        MD=MovieData.loadMatFile(ML.movieDataFile_{i});
        
        % Load the results of the detection and tracking
        outputDirDetect=[MD.outputDirectory_ filesep 'EB3' filesep detectionMethod];
        outputDirTrack=[outputDirDetect filesep 'plustipTrackerio'];       
        tmp=load([outputDirTrack filesep 'track' filesep 'trackNewFormat.mat']);
        tracks=tmp.tracks;
        
        outputDirDist=[MD.outputDirectory_ filesep 'EB3PoleRef' filesep poleDetectionMethod filesep detectionMethod]
        sphericalCoord=load([outputDirDist filesep 'sphericalCoord.mat']);
        dist=load([outputDirDist filesep 'dist.mat']);
        
        %Spherical radius
        [azimuth,elevation,time,trackId,poleId]=sphericalDistribution(tracks,sphericalCoord.azimuth,sphericalCoord.elevation,sphericalCoord.rho,dist.poleId,radius);
        cumulAzi{k,i}=azimuth;
        cumulElev{k,i}=elevation;
        cumulPoleId{k,i}=poleId; 
        % lifetime
        lifetimesHist{k,i}=hist([tracks.lifetime],lifetimeBins);
        meanLft(k,i)=nanmean([tracks.lifetime]);     
        
        % total displacement
        sx=MD.pixelSize_;sz=MD.pixelSizeZ_;
        lengths=arrayfun(@(t) sum(( [sx*t.x(end),sx*t.y(end),sz*t.z(end)]- ... 
                                    [sx*t.x(1),  sx*t.y(1),  sz*t.z(1)  ]).^2).^0.5 ,tracks);
        meanLength(k,i)=nanmean(lengths);     
        
        % speed
        st=MD.timeInterval_;
        speed=arrayfun(@(t) mean(sum(( [sx*t.x(1:end-1),sx*t.y(1:end-1),sz*t.z(1:end-1)]- ... 
                                 [sx*t.x(2:end),  sx*t.y(2:end),  sz*t.z(2:end)  ]).^2).^0.5/st) ,tracks);
        meanSpeed(k,i)=nanmean(speed);       
    end
end
% meanSpeed(2,1)=NaN;
% meanLength(2,1)=NaN;

%% plot results
handles=setupFigure(3,1,3,'Name',xpName);

%lifetime
colors={'r','b','g','k'};
for i=1:size(lifetimesHist,1)
    plot(handles(1),lifetimeBins,vertcat(lifetimesHist{i,:}),colors{i});
    hold on;
    xlim(handles(1),[min(lifetimeBins) max(lifetimeBins)]);
end
hold off
xlabel(handles(1),'lifetime (frame)')
ylabel(handles(1),'count')

v=[];names=[];
for i=1:size(lifetimesHist,1)
    h = findobj('Color',colors{i});
    v = [v h(1)];
end
legend(v,conditionName{:});
hold off

% mean Length
boxplot(handles(1),meanLft',conditionName);
ylabel(handles(1),'mean lifetime (s)');

% mean Length
boxplot(handles(2),meanLength',conditionName);
ylabel(handles(2),'mean distance (nm)');

% mean speed
boxplot(handles(3),meanSpeed',conditionName);
ylabel(handles(3),'growth (nm/s)');



%% Elevation vs Azimuth
%CellIdx={[1:],[1:5],[1:3]};
for MLIdx=1:length(aMovieListArray)
    figure('name',conditionName{MLIdx},'Position', [100, 100, 1024, 1200]);
    nCell=length(CellIdx{MLIdx});
    %handles=setupFigure(nCell,2,2*nCell,'Name',conditionName{MLIdx},'AxesWidth',4,'AxesHeight',4);
    for cIdx=CellIdx{MLIdx}
        %    title('Growth repartition at 4 microns from the poles')
        g=subplot(nCell,4,cIdx*4-3);
        p = get(g,'position');
        p(4) = p(4)*1.1; % Add 10 percent to height
        p(3) = p(3)*1.1; % Add 10 percent to height
        %p(1) = p(1)*0.1; % Add 10 percent to height
        %p(2) = p(2)*0.1; % Add 10 percent to height
        set(g, 'position', p);
        plotTipIdx=(cumulElev{MLIdx,cIdx}<0)&(cumulPoleId{MLIdx,cIdx}==1);
        polar(cumulAzi{MLIdx,cIdx}(plotTipIdx),abs(-pi/2-cumulElev{MLIdx,cIdx}(plotTipIdx)),'r+');
        title(['Cell ' num2str(cIdx) ' astral Pole 1']);
        
        g=subplot(nCell,4,cIdx*4-2);
        p = get(g,'position');
        p(4) = p(4)*1.1; % Add 10 percent to height
        p(3) = p(3)*1.1; % Add 10 percent to height
        %p(1) = p(1)*0.1; % Add 10 percent to height
        %p(2) = p(2)*0.1; % Add 10 percent to height
        set(g, 'position', p);
        plotTipIdx=(cumulElev{MLIdx,cIdx}>0)&(cumulPoleId{MLIdx,cIdx}==1);
        polar(cumulAzi{MLIdx,cIdx}(plotTipIdx),(pi/2-cumulElev{MLIdx,cIdx}(plotTipIdx)),'r+');
        title(['Cell ' num2str(cIdx) ' interpolar Pole 1']);       

        g=subplot(nCell,4,cIdx*4-1);
        p = get(g,'position');
        p(4) = p(4)*1.1; % Add 10 percent to height
        p(3) = p(3)*1.1; % Add 10 percent to height
        %p(1) = p(1)*0.1; % Add 10 percent to height
        %p(2) = p(2)*0.1; % Add 10 percent to height
        set(g, 'position', p);
        plotTipIdx=(cumulElev{MLIdx,cIdx}>0)&(cumulPoleId{MLIdx,cIdx}==2);
        polar(cumulAzi{MLIdx,cIdx}(plotTipIdx),(pi/2-cumulElev{MLIdx,cIdx}(plotTipIdx)),'r+');
        title(['Cell ' num2str(cIdx) ' interpolar Pole 2']);      
        
        g=subplot(nCell,4,cIdx*4);
        p = get(g,'position');
        p(4) = p(4)*1.1; % Add 10 percent to height
        p(3) = p(3)*1.1; % Add 10 percent to height
        %p(1) = p(1)*0.1; % Add 10 percent to height
        %p(2) = p(2)*0.1; % Add 10 percent to height
        set(g, 'position', p);
        plotTipIdx=(cumulElev{MLIdx,cIdx}<0)&(cumulPoleId{MLIdx,cIdx}==2);
        polar(cumulAzi{MLIdx,cIdx}(plotTipIdx),abs(-pi/2-cumulElev{MLIdx,cIdx}(plotTipIdx)),'r+');
        title(['Cell ' num2str(cIdx) ' astral Pole 2']);               
    end
end

%% Integrated Elevation vs Azimuth
figure('name','Integrated Elevation vs Azimuth','Position', [100, 100, 1024, 1200]);
colors={'r+','b+','g+','k+','y+','r.'};
for MLIdx=1:length(aMovieListArray)
    
    g=subplot(length(aMovieListArray),4,MLIdx*4-3);
    p = get(g,'position');
    p(4) = p(4)*1.1; % Add 10 percent to height
    p(3) = p(3)*1.1; % Add 10 percent to height
    set(g, 'position', p);
    nCell=length(aMovieListArray(MLIdx).movieDataFile_);    
    for cIdx=CellIdx{MLIdx}
        plotTipIdxAstralPole1=(cumulElev{MLIdx,cIdx}<0)&(cumulPoleId{MLIdx,cIdx}==1);
        polar(cumulAzi{MLIdx,cIdx}(plotTipIdxAstralPole1),abs(-pi/2-cumulElev{MLIdx,cIdx}(plotTipIdxAstralPole1)),colors{cIdx});
        hold on;        
    end
    title([ conditionName{MLIdx} ' astral Pole 1'])
    
    g=subplot(length(aMovieListArray),4,MLIdx*4-2);
    p = get(g,'position');
    p(4) = p(4)*1.1; % Add 10 percent to height
    p(3) = p(3)*1.1; % Add 10 percent to height
    set(g, 'position', p);
    for cIdx=CellIdx{MLIdx}
        plotTipIdxInterPole1=(cumulElev{MLIdx,cIdx}>0)&(cumulPoleId{MLIdx,cIdx}==1);
        polar(cumulAzi{MLIdx,cIdx}(plotTipIdxInterPole1),(pi/2-cumulElev{MLIdx,cIdx}(plotTipIdxInterPole1)),colors{cIdx});
        hold on;
    end
    title([ conditionName{MLIdx} ' interpolar Pole 1'])
    
    g=subplot(length(aMovieListArray),4,MLIdx*4-1);
    p = get(g,'position');
    p(4) = p(4)*1.1; % Add 10 percent to height
    p(3) = p(3)*1.1; % Add 10 percent to height
    set(g, 'position', p);
    for cIdx=CellIdx{MLIdx}
        plotTipIdxInterPole2=(cumulElev{MLIdx,cIdx}>0)&(cumulPoleId{MLIdx,cIdx}==2);

        polar(cumulAzi{MLIdx,cIdx}(plotTipIdxInterPole2)-pi,(pi/2-cumulElev{MLIdx,cIdx}(plotTipIdxInterPole2)),colors{cIdx});
        hold on;
    end
    title([ conditionName{MLIdx} ' interpolar Pole 2'])
    
    g=subplot(length(aMovieListArray),4,MLIdx*4);
    p = get(g,'position');
    p(4) = p(4)*1.1; % Add 10 percent to height
    p(3) = p(3)*1.1; % Add 10 percent to height
    set(g, 'position', p);
    nCell=length(aMovieListArray(MLIdx).movieDataFile_);    
    for cIdx=CellIdx{MLIdx}
        plotTipIdxAstralPole2=(cumulElev{MLIdx,cIdx}<0)&(cumulPoleId{MLIdx,cIdx}==2);
        polar((cumulAzi{MLIdx,cIdx}(plotTipIdxAstralPole2)),abs(-pi/2-cumulElev{MLIdx,cIdx}(plotTipIdxAstralPole2)),colors{cIdx});
        hold on;
    end
    title([ conditionName{MLIdx} ' astral Pole 2'])
    
end

