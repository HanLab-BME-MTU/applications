%%% Description
% This script shows how to pull tracking results and provides some basic
% plotting of the track output. 

%% %% USER INPUT 
xpName='prometaphase'; % For convenience, Figure window will bear this name.

%% LOADING MOVIES INFO
% MovieList paths: 

%MLPath='/work/gdanuser/proudot/project/EB3-3D-track/data-analysis/four-phases/'
%MLPath='/project/cellbiology/gdanuser/december/philippe/+tipTracking/chemical-inhibition/analysis/anaphase/control/';

% MovieList FileName (the combination of condition you want to compare). 
%movieListFileNames={'controlAnaphase.mat'};
%movieListFileNames={'prometaphaseCells.mat'};
%movieListFileNames={'controlMetaphase.mat','10nMMetaphase.mat'};
%movieListFileNames={'controlMetaphase.mat','10nMMetaphase.mat','33nMMetaphase.mat','100nMMetaphase.mat'};
%movieListFileNames={'controlInterphase.mat','10nMInterphase.mat','33nMInterphase.mat','100nMInterphase.mat',};
%movieListFileNames={'controlAnaphase.mat','10nMAnaphase.mat','33nMAnaphase.mat'};
MLPath='/project/bioinformatics/Danuser_lab/externBetzig/analysis/proudot/anaProject/sphericalProjection/prometaphase/cell1_12_half_volume_double_time/';
movieListFileNames={'movieList.mat'};
% Name associated to each MovieList
conditionName={'prometaphase'};
%conditionName={'control','10nM'};
%conditionName={'control','10nM','33nM','100nM'};

% Build the array of MovieList (auto)
% Build the array of MovieList (automatic)
aMovieListArray=cellfun(@(x) MovieList.loadMatFile([MLPath filesep x]), movieListFileNames,'unif',0);
aMovieListArray=aMovieListArray{:};


% Get the Cell index per condition
CellIdx=cell(1,length(movieListFileNames));
for i=1:length(CellIdx)
    CellIdx{i}=[1:length(aMovieListArray(i).movieDataFile_)];
end;


%% Input parameters
sphericalProjectionRadius=5000;
interpolarMTAngle=0.5;
detectionMethodEB3='pointSourceAutoSigmaLM';
detectionMethodKin='pointSourceAutoSigmaFit';  
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

% Azimuth and elevation for EB3
cumulAzi=cell(length(aMovieListArray),maxCellNb); 
cumulElev=cell(length(aMovieListArray),maxCellNb); 
cumulPoleId=cell(length(aMovieListArray),maxCellNb); 
cumulTimePt=cell(length(aMovieListArray),maxCellNb); 

% Azimuth and elevation for EB3 that caught a Kinetochore. 
EB3CatchingId=cell(length(aMovieListArray),maxCellNb); 

% Azimuth and elevation for Kin
cumulAziKin=cell(length(aMovieListArray),maxCellNb); 
cumulRhoKin=cell(length(aMovieListArray),maxCellNb); 
cumulTrackIdKin=cell(length(aMovieListArray),maxCellNb); 
cumulElevKin=cell(length(aMovieListArray),maxCellNb); 
cumulTimePtKin=cell(length(aMovieListArray),maxCellNb); 


%% Loop over each different conditions to collect all the data. 
for k=1:length(aMovieListArray)
    ML=aMovieListArray(k);
    %[~,conditionName{k}]=fileparts(fileparts(ML.movieListPath_));
    %% loop over each different cell in each condition
    for i=1:length(ML.movieDataFile_)
        MD=MovieData.loadMatFile(ML.movieDataFile_{i});
        
        % Load the tracking results 
        outputDirTrack=[MD.outputDirectory_ filesep 'EB3' filesep 'trackSpindleRef' filesep ];     
        tmp=load([outputDirTrack filesep 'tracksStageRef.mat']);
        tracks=tmp.tracksStageRef;
        
        % Load the spherical intersection
        outpurDir=[MD.outputDirectory_ filesep 'EB3' filesep 'sphericalProjection' filesep 'radius-' num2str(sphericalProjectionRadius)];
        sphericalProjection=load([outpurDir filesep 'sphericalProjection.mat']);
        cumulAzi{k,i}=sphericalProjection.sphericalAzimuth;
        cumulElev{k,i}=sphericalProjection.sphericalElevation;
        cumulPoleId{k,i}=sphericalProjection.poleId; 
        cumulTimePt{k,i}=sphericalProjection.time; 

        % Load the spherical coordinate of detected Kinetochore and their associated
        % track ID. 
        outputDirDetect=[MD.outputDirectory_ filesep 'Kin'  filesep 'detection' filesep];
        if(exist([outputDirDetect filesep 'sphericalCoordBothPoles.mat'], 'file') == 2)
            sphericalCoord=load([outputDirDetect filesep 'sphericalCoordBothPoles.mat']);

            outputDirTrack=[MD.outputDirectory_ filesep 'Kin' filesep 'track' filesep ];
            trackData=load([outputDirTrack  filesep 'tracksStageRef.mat']);
            
              sphCoordCumulKin=sphericalRadiusBinning(sphericalCoord.sphCoord,[0,sphericalProjectionRadius,100000],MD.timeInterval_,trackData.tracksStageRef,[]);
              cumulAziKin{k,i}=sphCoordCumulKin{2}.azimuth;
              cumulElevKin{k,i}=sphCoordCumulKin{2}.elevation;
              cumulRhoKin{k,i}=sphCoordCumulKin{2}.rho;
              cumulTimePtKin{k,i}=sphCoordCumulKin{2}.time;
              cumulTrackIdKin{k,i}=sphCoordCumulKin{2}.trackId;
            %%
            
            for fIdx=2:(length(sphericalCoord.sphCoord.elevation)-1)
                %%
                crossingAtFrameIdx=ceil(sphericalProjection.time)==fIdx;
                EB3Pos=[  sphericalProjection.sphericalAzimuth(crossingAtFrameIdx)' ...
                          sphericalProjection.sphericalElevation(crossingAtFrameIdx)'];
                EB3PoleId=sphericalProjection.poleId(crossingAtFrameIdx);
                EB3PosPole1=EB3Pos(EB3PoleId==1,:);
                EB3PosPole2=EB3Pos(EB3PoleId==2,:);
                FarKinIdx=(min(sphericalCoord.sphCoord.rho{fIdx},[],2)>sphericalProjectionRadius);
                KinPosPole1=[sphericalCoord.sphCoord.azimuth{fIdx}(FarKinIdx,1) sphericalCoord.sphCoord.elevation{fIdx}(FarKinIdx,1)];
                KinPosPole2=[sphericalCoord.sphCoord.azimuth{fIdx}(FarKinIdx,2) sphericalCoord.sphCoord.elevation{fIdx}(FarKinIdx,2)];
                [catchingEB3Idx1,caughtKinIdx1]=colocalizationLAP(EB3PosPole1,KinPosPole1,0.2);
                [catchingEB3Idx2,caughtKinIdx2]=colocalizationLAP(EB3PosPole2,KinPosPole2,0.2);
                %%
                EB3FrameIdxP1=find(EB3PoleId==1);EB3FrameIdxP1=EB3FrameIdxP1(catchingEB3Idx1);
                EB3FrameIdxP2=find(EB3PoleId==2);EB3FrameIdxP2=EB3FrameIdxP2(catchingEB3Idx2);
                EB3FrameIdx=[EB3FrameIdxP1 EB3FrameIdxP2];
                EB3Idx=find(crossingAtFrameIdx); EB3Idx=EB3Idx(EB3FrameIdx);
                EB3CatchingId{k,i}=[EB3CatchingId{k,i} EB3Idx ];
            end

              %%  
                   
            

        end
        
        % lifetime
        lifetimesHist{k,i}=hist([tracks.lifetime],lifetimeBins);
        meanLft(k,i)=nanmean([tracks.lifetime]);     
        
        % interpolar lifetime 
        interpolarTrackId=sphericalProjection.trackId(sphericalProjection.sphericalElevation<interpolarMTAngle);
        meanInterpolarLft(k,i)=nanmean([tracks(interpolarTrackId).lifetime]);     
        
        
        % total displacement
        sx=MD.pixelSize_;sz=MD.pixelSizeZ_;
        lengths=arrayfun(@(t) sum(( [sx*t.x(end),sx*t.y(end),sz*t.z(end)]- ... 
                                    [sx*t.x(1),  sx*t.y(1),  sz*t.z(1)  ]).^2).^0.5 ,tracks);
        meanLength(k,i)=nanmean(lengths);     
        
        % speed
        st=MD.timeInterval_;
        speed=arrayfun(@(t) nanmedian(nansum(([sx*t.x(1:end-1);sx*t.y(1:end-1);sz*t.z(1:end-1)]-[sx*t.x(2:end);  sx*t.y(2:end);  sz*t.z(2:end)]).^2).^0.5),tracks);

%         speed=arrayfun(@(t) mean(sum(( [sx*t.x(1:end-1),sx*t.y(1:end-1),sz*t.z(1:end-1)]- ... 
%                                  [sx*t.x(2:end),  sx*t.y(2:end),  sz*t.z(2:end)  ]).^2).^0.5/st) ,tracks);
        meanSpeed(k,i)=nanmean(speed);       
    end
end
% meanSpeed(2,1)=NaN;
% 'filled'
% meanLength(2,1)=NaN;

%% plot results
handles=setupFigure(3,1,3,'Name',xpName);

%lifetime
% colors={'r','b','g','k'};
% for i=1:size(lifetimesHist,1)
%     plot(handles(1),lifetimeBins,vertcat(lifetimesHist{i,:}),colors{i});
%     hold on;
%     xlim(handles(1),[min(lifetimeBins) max(lifetimeBins)]);
% end
% hold off
% xlabel(handles(1),'lifetime (frame)')
% ylabel(handles(1),'count')
% v=[];names=[];
% for i=1:size(lifetimesHist,1)
%     h = findobj('Color',colors{i});
%     v = [v h(1)];
% end
% legend(v,conditionName{:});
% hold off


% mean Length
boxplot(handles(1),meanLft',conditionName);
ylabel(handles(1),'mean lifetime (s)');


% mean Length
boxplot(handles(2),meanLength',conditionName);
ylabel(handles(2),'mean distance (nm)');

% mean speed
boxplot(handles(3),meanSpeed',conditionName);
ylabel(handles(3),'growth (nm/s)');


% mean Liftime interpolar
boxplot(handles(1),meanLft',conditionName);
ylabel(handles(1),'mean lifetime (s)');

%% plot results intepolar only
handles=setupFigure(3,1,3,'Name',xpName);
% mean Length
boxplot(handles(1),meanInterpolarLft',conditionName);
ylabel(handles(1),'mean lifetime (s)');

%% Elevation vs Azimuth
EB3MarkerSize=10;
KinMarkerSize=50;
cmapKin=jet(600);
cmapEB3=summer(150); %cmapEB3=cmapEB3(1:40,:);
temporalWindow=1; %Number of Frames used for integration
%CellIdx={[1:],[1:5],[1:3]};
for MLIdx=1:length(aMovieListArray)
    ML=aMovieListArray(k); 
    for cIdx=CellIdx{MLIdx}
        MD=MovieData.loadMatFile(ML.movieDataFile_{cIdx});
        handles=setupFigure(1,4,'Name',[ conditionName{MLIdx} ' Cell ' num2str(cIdx) ],'AxesWidth',4,'AxesHeight',4,'DisplayMode', 'print');
        plotSpindleSphericalProjection(handles,cumulAzi{MLIdx,cIdx},cumulElev{MLIdx,cIdx},cumulPoleId{MLIdx,cIdx},cumulTimePt{MLIdx,cIdx}, ...
            cumulAziKin{MLIdx,cIdx},cumulElevKin{MLIdx,cIdx},cumulTimePtKin{MLIdx,cIdx},cumulTrackIdKin{MLIdx,cIdx},[0 MD.nFrames_*MD.timeInterval_],EB3MarkerSize,KinMarkerSize,cmapEB3,cmapKin,[])
        outpurDir=[MD.outputDirectory_ filesep 'EB3' filesep 'sphericalProjection' filesep 'radius-' num2str(sphericalProjectionRadius)];
        mkdir(outpurDir);
        print([outpurDir filesep 'cumulative.eps'],'-depsc');
        mkdir([outpurDir filesep 'png']);
        mkdir([outpurDir filesep 'eps']);
        for t=1:temporalWindow:(MD.nFrames_-temporalWindow)
            [handles,~,fhandle]=setupFigure(1,4,'Name',[ conditionName{MLIdx} ' Cell ' num2str(cIdx) ],'AxesWidth',4,'AxesHeight',4,'DisplayMode', 'print');
            set(fhandle,'Visible','off');
            plotSpindleSphericalProjection(handles,cumulAzi{MLIdx,cIdx},cumulElev{MLIdx,cIdx},cumulPoleId{MLIdx,cIdx},cumulTimePt{MLIdx,cIdx}, ... 
                cumulAziKin{MLIdx,cIdx},cumulElevKin{MLIdx,cIdx},cumulTimePtKin{MLIdx,cIdx},cumulTrackIdKin{MLIdx,cIdx}, ... 
                MD.timeInterval_*[t t+temporalWindow],EB3MarkerSize,KinMarkerSize,cmapEB3,cmapKin, ...
                EB3CatchingId{MLIdx,cIdx});
            
            print([outpurDir filesep 'png' filesep 'time-' num2str(t,'%03d') '-' num2str(t+temporalWindow,'%03d') '.png'],'-dpng');
            print([outpurDir filesep 'eps' filesep 'time-' num2str(t,'%03d') '-' num2str(t+temporalWindow,'%03d') '.eps'],'-depsc');
        end 
    end
end

