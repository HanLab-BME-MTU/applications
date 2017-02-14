function printAllKinPoleRef(MD, kinTracks)
% $$$ ip = inputParser;
% $$$ ip.CaseSensitive = false;
% $$$ ip.KeepUnmatched = true;
% $$$ ip.addRequired('MD',@(MD) isa(MD,'MovieData'));
% $$$ ip.addParameter('kinTracks',[],@(x) isa(x,'Tracks'));
% $$$ ip.addParameter('printAll',false, @islogical);
% $$$ ip.parse(kinTrack,EB3Tracks,varargin{:});
% $$$ p=ip.Results;

outputDirProj=[MD.outputDirectory_ filesep 'Kin' filesep 'appearingMT' filesep 'projection'];
outputDirProjP1=[outputDirProj filesep 'P1' filesep];
outputDirProjP2=[outputDirProj filesep 'projection' filesep 'P2' filesep];
mkdir(outputDirProjP1);
mkdir(outputDirProjP2);
for kinIdx=1:length(kinTracks)
    kinTrack=kinTracks(kinIdx);
    progressText(kinIdx/length(kinTracks),'Print KinPole ref.');
    [handles,fhandle]=printKinPoleRef(kinTrack.KP1,kinTrack.appearingMTKinP1Ref);
    print([outputDirProjP1 'kin' num2str(kinIdx,'%03d') '.png'],'-dpng');
    hold(handles(1),'off');
    close(fhandle);

    [handles,fhandle]=printKinPoleRef(kinTrack.KP2,kinTrack.appearingMTKinP2Ref);
    print([outputDirProjP2 'kin' num2str(kinIdx,'%03d') '.png'],'-dpng');
    hold(handles(1),'off');
    close(fhandle);
end

