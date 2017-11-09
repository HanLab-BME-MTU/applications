function printAllDetectKinPoleRef(kinTracks,varargin)
ip = inputParser;
ip.CaseSensitive = false;
ip.KeepUnmatched = true;
ip.addRequired('kinTracks',@(x) isa(x,'Tracks'));
ip.addOptional('process',[],@(MD) isa(MD,'Process'));
ip.addOptional('name',[]);
ip.addOptional('colorPropP1',[]);
ip.addOptional('colorPropP2',[]);
ip.addOptional('colorPropRange',[]);
ip.parse(kinTracks,varargin{:});
p=ip.Results;

if(~isempty(p.process))
    outputDirProj=[p.process.getOwner.outputDirectory_ filesep 'Kin' filesep 'MTTipsBias' filesep p.name filesep 'projection'];
    mkdir2016a(outputDirProj);
end

for kinIdx=1:length(kinTracks)
    kinTrack=kinTracks(kinIdx);
    progressText(kinIdx/length(kinTracks),'Print KinPole ref.');
    [handles,~,fhandle]=setupFigure(1,2,2,'AxesWidth',8,'AxesHeight',4);
    if(~isempty(p.colorPropP1))
        %[handles,fhandle]=printKinPoleRef(kinTrack.KP1,kinTrack.appearingMTKinP1Ref,'EB3ColorAfter',getfield(kinTrack,p.colorPropP1),'valueRange',p.colorPropRange);
    else
        printDetectPoleRef(kinTrack.KP1,kinTrack.associatedDetectP1KinRef,'handle',handles(1))
    end

    if(~isempty(p.colorPropP2))
        %[handles,fhandle]=printKinPoleRef(kinTrack.KP2,kinTrack.appearingMTKinP2Ref,'EB3ColorAfter',getfield(kinTrack,p.colorPropP2),'valueRange',p.colorPropRange);
    else
        printDetectPoleRef(kinTrack.KP2,kinTrack.associatedDetectP2KinRef,'handle',handles(2))
    end
    if(~isempty(p.process))
        print([outputDirProj filesep 'kin' num2str(kinIdx,'%03d') '.png'],'-dpng');
        hold(handles(1),'off');
        close(fhandle);
    end

end

function mkdir2016a(dir)
system(['mkdir -p ' dir]);

