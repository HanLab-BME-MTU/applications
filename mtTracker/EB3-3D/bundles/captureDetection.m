
function kinTracks=captureDetection(MD,varargin)
% EB3 and Kin tracks need to be augmented with spherical coordinate and set
% in nanometers (function addSpindleRef.m)

ip = inputParser;
ip.CaseSensitive = false;
ip.KeepUnmatched = true;
ip.addRequired('MD',@(MD) isa(MD,'MovieData'));
ip.addParameter('kinTracks',[],@(x) isa(x,'Tracks'));
ip.addParameter('EB3tracks',[],@(x) isa(x,'Tracks'));
ip.addParamValue('name','',@ischar);
ip.addParameter('printAll',false, @islogical);
ip.addParameter('process',[]);
ip.addParameter('testKinIdx',[1],@isnumeric);
%ip.addParameter('testKinIdx',[19 46 156],@isnumeric);
ip.addParameter('distanceCutOff',0.1,@isnumeric);
ip.parse(MD,varargin{:});
p=ip.Results;

%% Load EB3 tracks add azimuth info, change coordinate to real space measurement.
%%

if(isempty(p.kinTracks)||isempty(p.EB3tracks))
    outputDirCatchingMT=[MD.outputDirectory_ filesep 'Kin' filesep 'track'];
    tmp=load([outputDirCatchingMT filesep 'augmentedSpindleRef.mat'],'kinTracks');
    kinTracks=tmp.kinTracks;
    outputDirCatchingMT=[MD.outputDirectory_ filesep 'EB3' filesep 'track'];
    tmp=load([outputDirCatchingMT filesep 'augmentedSpindleRef.mat'],'EB3Tracks');
    EB3tracks=tmp.EB3Tracks;
else
    kinTracks=p.kinTracks;
    EB3tracks=p.EB3tracks;
end



%%  For each Kinetochore, indentify the +TIP that disappear close to the Pole-Kin axis.
% The pole ID use for the poleKin axis is defined by the location of the first point of each
% track

testKinIdx=p.testKinIdx;
distanceCutOff=p.distanceCutOff;

EB3TermFrames=[EB3tracks.endFrame];
EB3TermPoleId=arrayfun(@(t) t.poleId(end),EB3tracks);
EB3TermAzi=arrayfun(@(t,p) t.azimuth(p,end),EB3tracks,EB3TermPoleId);
EB3TermRho=arrayfun(@(t,p) t.rho(p,end),EB3tracks,EB3TermPoleId);
EB3TermElev=arrayfun(@(t,p) t.elevation(p,end),EB3tracks,EB3TermPoleId);
for kIdx=1:length(kinTracks)
    progressText(kIdx/length(kinTracks),'Catching MT.');
    kinTrack=kinTracks(kIdx);
    try
        kinTrack.addprop('catchingMT');
    catch
    end;
    kinTrack.catchingMT=[];
    for pIdx=1:length(kinTrack.f)
        fIdx=kinTrack.f(pIdx);
        coexistingEB3=(EB3TermFrames==fIdx); % Co existence is characterized if EB3 disappear when Kin is still alive.
        if(any(coexistingEB3))
            MTPoles=EB3TermPoleId(coexistingEB3);
            aziDiff=abs(kinTrack.azimuth(MTPoles,pIdx) - EB3TermAzi(coexistingEB3));
            elevDiff=(abs(kinTrack.elevation(MTPoles,pIdx)- EB3TermElev(coexistingEB3)));
%             kinEB3DistP1=sqrt(sum(([min(aziDiff,2*pi-aziDiff) ...
%                                         (kinTrack.elevation(MTPoles,pIdx)- EB3TermElev(coexistingEB3))]).^2,2));
%             caughtMTIndex=find(coexistingEB3);
%             caughtMTIndex=caughtMTIndex(kinEB3DistP1<distanceCutOff);
            kinEB3DistP1=(min(aziDiff,2*pi-aziDiff)<distanceCutOff) & (min(elevDiff,pi-elevDiff)<distanceCutOff);% ...
                          %& (EB3TermRho(coexistingEB3)<kinTrack.rho(MTPoles,pIdx)*1.02); % Rho is below kin row with a 2% margin. 
            caughtMTIndex=find(coexistingEB3);
            caughtMTIndex=caughtMTIndex(kinEB3DistP1);
            if(~isempty(caughtMTIndex))
                kinTrack.catchingMT=[kinTrack.catchingMT; EB3tracks(caughtMTIndex)] ;
                for mtIdx=caughtMTIndex
                    try
                        EB3tracks(mtIdx).addprop('caughtKin');
                    catch
                    end;
                    EB3tracks(mtIdx).caughtKin=[EB3tracks(mtIdx).caughtKin kinTrack];
                end;
            end

        end
    end
end

%% For test kinetochore, plot an Amira file with attached mt
outputDirAmira=[MD.outputDirectory_ filesep 'Kin' filesep 'catchingMT' filesep p.name filesep 'Test' filesep 'Amira'];
for kIdx=min(length(testKinIdx),testKinIdx)
    kinTrack=kinTracks(kIdx);
    trackSet=[kinTrack; kinTrack.catchingMT];
    trackType=[1; zeros(length(kinTrack.catchingMT),1)];
    amiraWriteTracks([outputDirAmira filesep 'kin_' num2str(kIdx) filesep 'kin_' num2str(kIdx)  '.am'],trackSet,'edgeProp',{{'kinEB3',trackType}})
end
%%
if(p.printAll)
    %% For each kinetochore, plot an Amira file with attached mt
    outputDirAmira=[MD.outputDirectory_ filesep 'Kin' filesep 'catchingMT' filesep p.name filesep  'Amira' filesep];
    for kIdx=1:length(kinTracks)
        kinTrack=kinTracks(kIdx);
        trackSet=[kinTrack; kinTrack.catchingMT];
        trackType=[1; zeros(length(kinTrack.catchingMT),1)];
        dataIsotropy=[MD.pixelSize_ MD.pixelSize_ MD.pixelSize_];
        amiraWriteTracks([outputDirAmira filesep 'kin_' num2str(kIdx) '.am'],trackSet,'cumulativeOnly',true,'edgeProp',{{'kinEB3',trackType}})
    end
end


%% Load the pole info
poleDetectionMethod=['simplex_scale_' num2str(3,'%03d')];
outputDirPoleDetect=[MD.outputDirectory_ filesep 'EB3' filesep 'poles' filesep poleDetectionMethod filesep];
tmp=load([outputDirPoleDetect filesep 'poleDetection.mat']);
poleMovieInfo=tmp.poleMovieInfo;

% WARNING: this is not a trajectory, merely a collection of poles to ease
% implementation.
P1=struct();
P1.x=arrayfun(@(d) MD.pixelSize_*(d.xCoord(1,1)-1)+1,poleMovieInfo)';
P1.y=arrayfun(@(d) MD.pixelSize_*(d.yCoord(1,1)-1)+1,poleMovieInfo)';
P1.z=arrayfun(@(d) MD.pixelSize_*(d.zCoord(1,1)-1)+1,poleMovieInfo)';

P2=struct();
P2.x=arrayfun(@(d) MD.pixelSize_*(d.xCoord(2,1)-1)+1,poleMovieInfo)';
P2.y=arrayfun(@(d) MD.pixelSize_*(d.yCoord(2,1)-1)+1,poleMovieInfo)';
P2.z=arrayfun(@(d) MD.pixelSize_*(d.zCoord(2,1)-1)+1,poleMovieInfo)';

poleStructs=[P1 P2];

%% For each MT captured by a kinetochore, create a MT in the referential defined by kin Pole referential.
% the z axis is the kinetochore-pole axis
% the x axis is
for kIdx=1:length(kinTracks)
    progressText(kIdx/length(kinTracks),'2D projection.');
    %%
    kinTrack=kinTracks(kIdx);
    % Vector basis w.r.t each poles and kinetochore

    KP1=struct();
    KP1.x=kinTrack.x-P1.x(kinTrack.f);
    KP1.y=kinTrack.y-P1.y(kinTrack.f);
    KP1.z=kinTrack.z-P1.z(kinTrack.f);
    KP1.vZ=[KP1.x; KP1.y; KP1.z];
    KP1.vZ=KP1.vZ./repmat(sum(KP1.vZ.^2,1).^0.5,3,1);
    KP1.vX=[0*KP1.vZ(1,:);KP1.vZ(3,:);-KP1.vZ(2,:)];
    KP1.vY=cross(KP1.vX,KP1.vZ,1);

    KP2=struct();
    KP2.x=kinTrack.x-P2.x(kinTrack.f);
    KP2.y=kinTrack.y-P2.y(kinTrack.f);
    KP2.z=kinTrack.z-P2.z(kinTrack.f);
    KP2.vZ=[KP2.x; KP2.y; KP2.z];
    KP2.vZ=KP2.vZ./repmat(sum(KP2.vZ.^2,1).^0.5,3,1);
    KP2.vX=[0*KP2.vZ(1,:);KP2.vZ(3,:);-KP2.vZ(2,:)];
    KP2.vY=cross(KP2.vX,KP2.vZ,1);
    %%
    kinAxis=[KP1 KP2];

    % For each kinetochore, it would be elegant to create two version on the kinetochore-pole
    % However is the info is yield by the rho.

    % Register in original MT
    try
        kinTrack.addprop('catchingMTKinRef');
    catch
    end;
    kinTrack.catchingMTKinRef=kinTrack.catchingMT;
    for mIdx=1:length(kinTrack.catchingMT)
        mt=kinTrack.catchingMT(mIdx);
        basis=kinAxis(mt.poleId(1));

        % Copying EB3 track
        mtKinRef=mt.copy();
        mtKinRef.addprop('refname');
        mtKinRef.refname='Kin';
        mtKinRef.addprop('basis');
        mtKinRef.basis=basis;

        % Register in original MT
        try
            mt.addprop('altRef');
            kinTrack.addprop('kin');
        catch
        end;
        mt.altRef=[mt.altRef mtKinRef];
        kinTrack.catchingMTKinRef(mIdx)=mtKinRef;
        for pIdx=1:length(mt.f)
            f=min(mt.f(pIdx),length(basis.vX));
            B=[basis.vX(:,f) basis.vY(:,f) basis.vZ(:,f)];
            P=poleStructs(mt.poleId(1));
            recentered=[(mt.x(pIdx)-P.x(f)) (mt.y(pIdx)-P.y(f)) (mt.z(pIdx)-P.z(f))];
            v=recentered*B;
            mtKinRef.x(pIdx)=v(1); mtKinRef.y(pIdx)=v(2); mtKinRef.z(pIdx)= v(3);
        end;

    end
end

%%

process=ip.Results.process;
if(~isempty(process))
    outputDirCatchingMT=[MD.outputDirectory_ filesep 'Kin' filesep 'catchingMT' filesep p.name];
    mkdir(outputDirCatchingMT);
    save([outputDirCatchingMT filesep 'catchingMT.mat'],'kinTracks','EB3tracks')
    process.setOutFilePaths({[outputDirCatchingMT filesep 'catchingMT.mat']})    
    pa = process.getParameters();
    pa.parameters = p;
    process.setParameters(pa);
    process.setDateTime();
end 

%% First test, display the +TIP coordinate on a lateral view of the poleKin axis.

if(p.printAll)

end
