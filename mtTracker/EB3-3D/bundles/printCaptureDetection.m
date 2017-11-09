
function printCaptureDetection(MD,varargin)
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

outputDirProj=[MD.outputDirectory_ filesep 'Kin' filesep 'projections' filesep p.name filesep  'firstTest' filesep]
system(['mkdir ' outputDirProj]);

for kIdx=1:length(kinTracks)
%%
    kinTrack=kinTracks(kIdx);

    [handles,~,fhandle]=setupFigure(1,2,'AxesWidth',8,'AxesHeight',4,'DisplayMode', 'print');
    hold(handles(1),'on');
    hold(handles(2),'on');

    for poleId=1:2
        scatter(handles(poleId),kinTrack.rho(poleId,:),zeros(size(kinTrack.rho(poleId,:))),'r');
        scatter(handles(poleId),0,0,'g');
    end

    for mIdx=1:length(kinTrack.catchingMT)
        mt=kinTrack.catchingMT(mIdx);
        mtKinRef=kinTrack.catchingMTKinRef(mIdx);

        %Project on the plan defined by the poleKin axis and the interpolar
        %axis.
        plot(handles(mt.poleId(1)),mtKinRef.z,mtKinRef.x,'b-');
        ylim(handles(mt.poleId(1)),[-2000 2000])
        xlabel(handles(mt.poleId(1)),'Pole-Kinetochore axis (nm)')
        ylabel(handles(mt.poleId(1)),'Normal plane (nm)')
    end
%%


    print([outputDirProj 'kin' num2str(kIdx,'%03d') '.png'],'-dpng');

    hold(handles(1),'off');
    hold(handles(2),'off');
    close(fhandle);
end
