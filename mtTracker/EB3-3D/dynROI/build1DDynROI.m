function [dynROIs]=build1DDynROI(originsTracksOrProcess,ZTracksOrProcess,mappingDistance,varargin)
ip = inputParser;
ip.CaseSensitive = false;
ip.KeepUnmatched=true;
ip.addRequired('originsTracksOrProcess');
ip.addRequired('ZTracksOrProcess');
ip.parse(originsTracksOrProcess,ZTracksOrProcess,varargin{:});
p=ip.Results;

[refs,ROIs]=buildRefsFromTracks(originsTracksOrProcess,ZTracksOrProcess,varargin{:});
rIdx=1;
dynROI.ref=refs(rIdx);
dynROI.ROI=ROIs{rIdx};
dynROI.mappingDistance=mappingDistance;
dynROI.name='';
dynROIs(numel(refs))=dynROI;
for oIdx=1:size(refs,1)
	for zIdx=1:size(refs,2)
		dynROI.ref=refs(rIdx);
		dynROI.ROI=ROIs{rIdx};
		dynROI.mappingDistance=mappingDistance;
		dynROI.name=['O-' num2str(ROIs{rIdx}(1).index()) '-Z-' num2str(ROIs{rIdx}(2).index())];
		dynROIs(rIdx)=dynROI;
		rIdx=rIdx+1;
	end
end
