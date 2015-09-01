%%% Description
% This script shows how to pull tracking results and provides some basic
% plotting of the track output. 

%% %% USER INPUT 

%% INPUT MOVIES
% After movie loading (see loadMoviesManagementFile.m)
aMovieListArray=[MLControlAna MLnM10Ana MLnM33Ana];
conditionName={'control','10nM','33nM'}

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

% Per-cell avg track length 
meanLength=nan(length(aMovieListArray),maxCellNb);

% Per-cell avg track speed
meanSpeed=nan(length(aMovieListArray),maxCellNb);

% Azimuth and elevation
cumulAzi=cell(length(aMovieListArray),maxCellNb); 
cumulElev=cell(length(aMovieListArray),maxCellNb); 



%conditionName=cell(1,length(aMovieListArray));
%% loop over each different conditions
for k=1:length(aMovieListArray)
    ML=aMovieListArray(k);
    %[~,conditionName{k}]=fileparts(fileparts(ML.movieListPath_));
    %% loop over each different cell in each condition
    for i=1:length(ML.movieDataFile_)
        MD=ML.getMovie(i);
        
        % Load the results of the detection and tracking
        outputDirDetect=[MD.outputDirectory_ filesep 'EB3' filesep detectionMethod];
        outputDirTrack=[outputDirDetect filesep 'plustipTrackerio'];       
        tmp=load([outputDirTrack filesep 'track' filesep 'trackNewFormat.mat']);
        tracks=tmp.tracks;
        outputDirDist=[MD.outputDirectory_ filesep 'EB3PoleRef' filesep poleDetectionMethod filesep detectionMethod]
        sphericalCoord=load([outputDirDist filesep 'sphericalCoord.mat']);
        
        % Spherical radius
%         [azimuth,elevation,time,trackId]=sphericalDistribution(tracks,sphericalCoord.azimuth,sphericalCoord.elevation,sphericalCoord.rho,radius);
%         cumulAzi{k,i}=azimuth;
%         cumulElev{k,i}=elevation;
         
        % lifetime
        lifetimesHist{k,i}=hist([tracks.lifetime],lifetimeBins);
        
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
meanSpeed(2,1)=NaN;
meanLength(2,1)=NaN;

%% plot results
handles=setupFigure(3,1,3);

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
boxplot(handles(2),meanLength',conditionName);
ylabel(handles(2),'mean length');

% mean speed
boxplot(handles(3),meanSpeed',conditionName);
ylabel(handles(3),'growth');

%% Elevation vs Azimuth
figure();
title('Growth repartition at 4 microns from the poles')
subplot(1,2,1);
polar(cumulAzi{1,2}(cumulElev{1,2}<0),abs(-pi/2-cumulElev{1,2}(cumulElev{1,2}<0)),'r+');
title('astral');

subplot(1,2,2);
polar(cumulAzi{1,2}(cumulElev{1,2}>0),pi/2-cumulElev{1,2}(cumulElev{1,2}>0),'r+');
title('interpolar');

