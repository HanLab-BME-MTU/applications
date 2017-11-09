function dynROIs=renderPKTROI(dynROIs,varargin)
ip = inputParser;
ip.CaseSensitive = false;
ip.KeepUnmatched=true;
ip.addRequired('process');
ip.parse(originsTracksOrProcess,ZTracksOrProcess,varargin{:});
p=ip.Results;

for dIdx=1:length(dynROIs)
    tIdx=projKin(dIdx);
    pIdx=projPole(dIdx);
    if(isempty(dynROIs(dIdx))
        tic;
        dynROI=dynROIs(dIdx);
        track=dynROI.ROI(2);
        PTrack=dynROI.ROI(1);
        processProj=ExternalProcess(MD,'dynROIProj');
        processRenderer=ExternalProcess(MD,'dynROIProj');
        ref=dynROI.ref;
        projDynROI(MD,dynROIs(dIdx).ROI,'dynPoligonREF',ref.applyBase(dynROIs(dIdx).ROI,''),'FoF',ref, ...
            'name',['PK-P' num2str(PTrack.index) '-' num2str(track.index], ...
            'channelRender','grayRed','processSingleProj',processProj,'processRenderer',processRenderer ...
            'intMinPrctil',[20 98],'intMaxPrctil',[100 100],'fringeWidth',30,'insetFringeWidth',dynROI.mappingDistance);
        dynROIs(dIdx).renderer=processProj;
        toc;
    end
end