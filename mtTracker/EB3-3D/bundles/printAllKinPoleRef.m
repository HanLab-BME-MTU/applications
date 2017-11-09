function printAllKinPoleRef(kinTracks,varargin)
ip = inputParser;
ip.CaseSensitive = false;
ip.KeepUnmatched = true;
ip.addRequired('kinTracks',@(x) isa(x,'Tracks'));
ip.addOptional('MD',[],@(MD) isa(MD,'MovieData'));
ip.addOptional('name',[]);
ip.addOptional('colorPropP1',[]);
ip.addOptional('colorPropP2',[]);
ip.addOptional('colorPropRange',[]);
ip.parse(kinTracks,varargin{:});
p=ip.Results;

if(~isempty(p.MD))
    outputDirProj=[p.MD.outputDirectory_ filesep 'Kin' filesep p.name filesep 'appearingMT' filesep 'projection'];
    outputDirProjP1=[outputDirProj filesep 'P1' filesep];
    outputDirProjP2=[outputDirProj filesep 'P2' filesep];
    mkdir2016a(outputDirProjP1);
    mkdir2016a(outputDirProjP2);
end
for kinIdx=1:length(kinTracks)
    kinTrack=kinTracks(kinIdx);
    progressText(kinIdx/length(kinTracks),'Print KinPole ref.');
    if(~isempty(p.colorPropP1))
        [handles,fhandle]=printKinPoleRef(kinTrack.KP1,kinTrack.appearingMTKinP1Ref,'EB3ColorAfter',getfield(kinTrack,p.colorPropP1),'valueRange',p.colorPropRange);
    else
        [handles,fhandle]=printKinPoleRef(kinTrack.KP1,kinTrack.appearingMTKinP1Ref);
    end

 
    if(~isempty(p.MD))
    print([outputDirProjP1 'kin' num2str(kinIdx,'%03d') '.png'],'-dpng');
    hold(handles(1),'off');
    close(fhandle);    
    end;
    
    if(~isempty(p.colorPropP2))
        [handles,fhandle]=printKinPoleRef(kinTrack.KP2,kinTrack.appearingMTKinP2Ref,'EB3ColorAfter',getfield(kinTrack,p.colorPropP2),'valueRange',p.colorPropRange);
    else
        [handles,fhandle]=printKinPoleRef(kinTrack.KP2,kinTrack.appearingMTKinP2Ref);
    end
    if(~isempty(p.MD))
        print([outputDirProjP2 'kin' num2str(kinIdx,'%03d') '.png'],'-dpng');
        hold(handles(1),'off');
        close(fhandle);
    end

end

function mkdir2016a(dir)
system(['mkdir -p ' dir]);

