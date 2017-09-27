function dynROIs=renderPKTDynROI(MD,dynROIs,varargin)
ip = inputParser;
ip.CaseSensitive = false;
ip.KeepUnmatched=true;
ip.addRequired('MD');
ip.addRequired('dynROIs');
ip.addParamValue('ref',[]);
ip.addParamValue('contextDynROI',[]);
ip.addParamValue('namePrefix',[]);
ip.parse(MD,dynROIs,varargin{:});
p=ip.Results;

for dIdx=1:length(dynROIs)
        tic;
        dynROI=dynROIs(dIdx);
        track=dynROI.ROI(2);
        PTrack=dynROI.ROI(1);
        name=dynROI.name;

        processProj=ExternalProcess(MD,'dynROIProj');
        processRenderer=ExternalProcess(MD,'dynROIProj');
        if(isempty(p.ref))
            ref=dynROI.ref;
        else
            ref=p.ref;
        end
        if(isempty(p.contextDynROI))
            context=ref.applyBase(dynROIs(dIdx).ROI,'');
            fringeWidth=30;
        else
            context=ref.applyBase(p.contextDynROI.ROI,'');
            fringeWidth=p.contextDynROI.mappingDistance;
        end
        projectDynROI(MD,dynROIs(dIdx).ROI,'dynPoligonREF',context,'FoF',ref, ...
            'name',[p.namePrefix name], ...
            'channelRender','grayRed','processSingleProj',processProj,'processRenderer',processRenderer, ...
            'intMinPrctil',[20 70],'intMaxPrctil',[99.99 99.99],'fringeWidth',fringeWidth,'insetFringeWidth',dynROI.mappingDistance);
        dynROIs(dIdx).renderer=processRenderer;
        toc;
end