MD=MovieData.loadMatFile('/project/bioinformatics/Danuser_lab/externBetzig/analysis/proudot/anaProject/sphericalProjection/prometaphase/cell1_12_half_volume_double_time/movieData.mat');

%% Load EB3 and collect EB3 depolymerization onset and associated track ID 
outputDirDetect=[MD.outputDirectory_ filesep 'EB3' filesep 'detection' filesep];mkdir(outputDirDetect);
EB3SphCoord=load([outputDirDetect filesep 'sphericalCoord.mat']);
EB3PoleId=load([outputDirDetect filesep 'dist.mat']);

outputDirTrack=[MD.outputDirectory_ filesep 'EB3' filesep 'track' filesep  ];
tmp=load([outputDirTrack filesep 'tracksSpindleRef.mat']);
EB3tracks=tmp.tracksSpindleRef;

trackEndIndex=cell(1,EB3tracks.numTimePoints);
trackEndTrackId=cell(1,EB3tracks.numTimePoints);


%% Augment the structures with spherical Coordinate. 

% For MT
for tIdx=1:length(EB3tracks)
    tr=EB3tracks(tIdx);
    try
        tr.addprop('azimuth');
        tr.addprop('elevation');
        tr.addprop('rho');
        tr.addprop('poleId');
    catch
    end;

    nonGap=~tr.gapMask();
    tr.azimuth=nan(size(tr.f));
    tr.azimuth(nonGap)=arrayfun(@(i,f) EB3SphCoord.sphCoordBest.azimuth{f}(i), tr.tracksFeatIndxCG(nonGap),tr.f(nonGap));
    tr.elevation=nan(size(tr.f));
    tr.elevation(nonGap)=arrayfun(@(i,f) EB3SphCoord.sphCoordBest.elevation{f}(i), tr.tracksFeatIndxCG(nonGap),tr.f(nonGap));
    tr.rho=nan(size(tr.f));
    tr.rho(nonGap)=arrayfun(@(i,f) EB3SphCoord.sphCoordBest.rho{f}(i), tr.tracksFeatIndxCG(nonGap),tr.f(nonGap));
    tr.poleId=nan(size(tr.f));
    tr.poleId(nonGap)=arrayfun(@(i,f) EB3PoleId.poleId{f}(i), tr.tracksFeatIndxCG(nonGap),tr.f(nonGap));
end

%% load Kinetochores spherical coordinate
outputDirDetect=[MD.outputDirectory_ filesep 'Kin'  filesep 'detection' filesep];
kinSphericalCoord=load([outputDirDetect filesep 'sphericalCoordBothPoles.mat']);
kinSphericalCoord=kinSphericalCoord.sphCoord;

outputDirTrack=[MD.outputDirectory_ filesep 'Kin' filesep 'track' filesep ];
kinTrackData=load([outputDirTrack  filesep 'tracksSpindleRef.mat']);
kinTracks=kinTrackData.tracksSpindleRef;

% For kin set the spindle ref as the 1 point referential.
poleID=1;
for tIdx=1:length(kinTracks)
    tr=kinTracks(tIdx);
    try
        tr.addprop('azimuth');
        tr.addprop('elevation');
        tr.addprop('rho');
    catch
    end;

    nonGap=~tr.gapMask();
    tr.azimuth=nan(size(tr.f));
    tr.azimuth(nonGap)=arrayfun(@(i,f) kinSphericalCoord.azimuth{f}(i,poleID), tr.tracksFeatIndxCG(nonGap),tr.f(nonGap));
    tr.elevation=nan(size(tr.f));
    tr.elevation(nonGap)=arrayfun(@(i,f) kinSphericalCoord.elevation{f}(i,poleID), tr.tracksFeatIndxCG(nonGap),tr.f(nonGap));
    tr.rho=nan(size(tr.f));
    tr.rho(nonGap)=arrayfun(@(i,f) kinSphericalCoord.rho{f}(i,poleID), tr.tracksFeatIndxCG(nonGap),tr.f(nonGap));
end


%% For each Kinetochore, indentify the +TIP that disappear close to the Pole-Kin axis.

distanceCutOff=0.2;

EB3TermFrames=[EB3tracks.endFrame];
EB3TermAzi=arrayfun(@(t) t.azimuth(end),EB3tracks);
EB3TermRho=arrayfun(@(t) t.rho(end),EB3tracks);
EB3TermElev=arrayfun(@(t) t.elevation(end),EB3tracks);
for kIdx=1:length(kinTracks)
    kinTrack=kinTracks(kIdx);
    try
        kinTrack.addprop('catchingMT');
    catch
    end;
    kinTrack.caughtEB3=[];
    for pIdx=1:length(kinTrack.f)
        fIdx=kinTrack.f(pIdx);
        kinDetIdx=kinTrack.tracksFeatIndxCG(pIdx);
        kinCoordP1=[kinSphericalCoord.azimuth{fIdx}(kinDetIdx,1)' kinSphericalCoord.elevation{fIdx}(kinDetIdx,1)'];
        if(any(EB3TermFrames==fIdx))
            kinEB3DistP1=norm([(kinCoordP1(1) - EB3TermAzi(EB3TermFrames==fIdx))' (kinCoordP1(2)- EB3TermElev(EB3TermFrames==fIdx))']);
            caughtMTIndex=find(EB3TermFrames==fIdx);
            caughtMTIndex=caughtMTIndex(kinEB3DistP1<distanceCutOff);
            if(~isempty(caughtMTIndex))
                kinTrack.catchingMT=[kinTrack.catchingMT EB3tracks(caughtMTIndex)] ;             
                try
                    EB3tracks(caughtMTIndex).addprop('caughtKin');
                catch
                end;
                EB3tracks(caughtMTIndex).caughtKin=[EB3tracks(caughtMTIndex).caughtKin deal(kinTrack)];
            end
            
        end
    end
end

%% For each MT captured by a kinetochore, create a MT in the referential defined by the MT end and the 
for kIdx=1:length(kinTracks)
    for mIdx=1:length(kinTrack.catchingMT)
        % Copying EB3 track and linking it to the original.
        mt=kinTracks.catchingMT(mIdx));
        mtKinRef=mt.copy();
        mtKinRef.addprop('refname'); 
        mtKinRef.refname='Kin';
       
        mt.addprop('altRef');
        mt.atlRef=[mt.atlRef mtKinRef];

        mtKinRef.azimuth=mtKinRef.azimuth-kinTracks.azimuth(end);
        mtKinRef.elevation=mtKinRef.elevation-kinTracks.elevation(end);     
    end
end

%% First test, display the +TIP coordinate on a lateral view of the poleKin axis. 
for kIdx=length(kinSphericalCoord.sphCoord.azimuth{fIdx})
    [handles,~,fhandle]=setupFigure(1,2,'AxesWidth',4,'AxesHeight',4,'DisplayMode', 'print');
    for mIdx=1:length(kinTrack.catchingMT)
        mt=kinTracks.catchingMT(mIdx));
        mt.altRef(1);
        plot(handles(1),)
        
        
        plotSphericalProjection(handles(1),azimuths,elevations,colorIndex,markerSize,cmap, varargin)
        print([outputDir 'time-' num2str(t,'%03d') '-' num2str(t+temporalWindow,'%03d') '.png'],'-dpng');
    end
end

            
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