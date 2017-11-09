sphericalProjectionRadius=5000;
MD=MovieData.loadMatFile('/project/bioinformatics/Danuser_lab/externBetzig/analysis/proudot/anaProject/sphericalProjection/prometaphase/cell1_12_half_volume_double_time/movieData.mat');

%% Load EB3 and collect EB3 depolymerization onset
outputDirDetect=[MD.outputDirectory_ filesep 'EB3' filesep 'detection' filesep];mkdir(outputDirDetect);
EB3SphCoord=load([outputDirDetect filesep 'sphericalCoord.mat']);
EB3PoleId=load([outputDirDetect filesep 'dist.mat']);

outputDirTrack=[MD.outputDirectory_ filesep 'EB3' filesep 'track' filesep  ];
tmp=load([outputDirTrack filesep 'tracksStageRef.mat']);
tracks=tmp.tracksStageRef;

trackEndIndex=cell(1,tracks.numTimePoints);
for tIdx=1:length(tracks)
    if(tracks(tIdx).lifetime>3)
    trackEndIndex{tracks(tIdx).endFrame}=[trackEndIndex{tracks(tIdx).endFrame} tracks(tIdx).tracksFeatIndxCG(end)];
    end
end

EB3SphCoordEnd.azimuth=cellfun(@(a,t) a(t),EB3SphCoord.sphCoordBest.azimuth,trackEndIndex,'unif',0);
EB3SphCoordEnd.elevation=cellfun(@(a,t) a(t),EB3SphCoord.sphCoordBest.elevation,trackEndIndex,'unif',0);
EB3SphCoordEnd.rho=cellfun(@(a,t) a(t),EB3SphCoord.sphCoordBest.rho,trackEndIndex,'unif',0        );
EB3poleIdEnd=cellfun(@(a,t) a(t),EB3PoleId.poleId,trackEndIndex,'unif',0        );


%% load Kinetochores spherical coordinate
outputDirDetect=[MD.outputDirectory_ filesep 'Kin'  filesep 'detection' filesep];
sphericalCoord=load([outputDirDetect filesep 'sphericalCoordBothPoles.mat']);

outputDirTrack=[MD.outputDirectory_ filesep 'Kin' filesep 'track' filesep ]
trackData=load([outputDirTrack  filesep 'tracksStageRef.mat']);

            
%% Bin each set of coordinate in a sphere.             
radii=[1000 4000 7000 9000];
%radius=[0,sphericalProjectionRadius,100000];
KinSphCoordCumulKin=sphericalRadiusBinning(sphericalCoord.sphCoord,radii,MD.timeInterval_,trackData.tracksStageRef,[]);
EB3SphCoordEndCumul=sphericalRadiusBinning(EB3SphCoordEnd,radii,MD.timeInterval_,[],EB3poleIdEnd);

%% display binning results
%% cumulative

EB3MarkerSize=10;
KinMarkerSize=50;
cmapKin=jet(600);
cmapEB3=summer(150); %cmapEB3=cmapEB3(1:40,:);
temporalWindow=1; %Number of Frames used for integration

handles=setupFigure(length(radii)-1,4,4*(length(radii)-1),'Name','test binning','AxesWidth',4,'AxesHeight',4,'DisplayMode', 'print');
for rIdx=1:(length(radii)-1)
    plotSpindleSphericalProjection(handles(((rIdx-1)*4+1):rIdx*4),EB3SphCoordEndCumul{rIdx}.azimuth,EB3SphCoordEndCumul{rIdx}.elevation,EB3SphCoordEndCumul{rIdx}.poleId,EB3SphCoordEndCumul{rIdx}.time, ...
        KinSphCoordCumulKin{rIdx}.azimuth,KinSphCoordCumulKin{rIdx}.elevation,KinSphCoordCumulKin{rIdx}.time,KinSphCoordCumulKin{rIdx}.trackId,[0 MD.nFrames_*MD.timeInterval_],EB3MarkerSize,KinMarkerSize,cmapEB3,cmapKin,[])
end

%% Time split
temporalWindow=1;
outputDir=[MD.outputDirectory_ filesep 'EB3' filesep 'sphericalProjection' filesep 'radii'  num2str(radii,'-%d') filesep num2str(temporalWindow) filesep];
mkdir(outputDir)

for t=1:temporalWindow:(MD.nFrames_-temporalWindow)/2
    [handles,~,fhandle]=setupFigure(length(radii)-1,4,'AxesWidth',4,'AxesHeight',4,'DisplayMode', 'print');
    set(fhandle,'Visible','off');
    for rIdx=1:(length(radii)-1)
        plotSpindleSphericalProjection(handles(((rIdx-1)*4+1):rIdx*4),EB3SphCoordEndCumul{rIdx}.azimuth,EB3SphCoordEndCumul{rIdx}.elevation,EB3SphCoordEndCumul{rIdx}.poleId,EB3SphCoordEndCumul{rIdx}.time, ...
            KinSphCoordCumulKin{rIdx}.azimuth,KinSphCoordCumulKin{rIdx}.elevation,KinSphCoordCumulKin{rIdx}.time,KinSphCoordCumulKin{rIdx}.trackId,MD.timeInterval_*[t t+temporalWindow],EB3MarkerSize,KinMarkerSize,cmapEB3,cmapKin,[])
    end
    print([outputDir 'time-' num2str(t,'%03d') '-' num2str(t+temporalWindow,'%03d') '.png'],'-dpng');
    %print([outputDir filesep 'png' filesep 'time-' num2str(t,'%03d') '-' num2str(t+temporalWindow,'%03d') '.png'],'-dpng');
    %print([outputDir filesep 'eps' filesep 'time-' num2str(t,'%03d') '-' num2str(t+temporalWindow,'%03d') '.eps'],'-depsc');
end