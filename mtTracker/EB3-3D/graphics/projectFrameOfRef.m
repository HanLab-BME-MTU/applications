function processList=projectFrameOfRef(MD,ref,varargin)
% dynPoligonREF define the ROI in the final ref and is used for the actual displayed ROI location.
ip = inputParser;
ip.CaseSensitive = false;
ip.KeepUnmatched = true;
ip.addRequired('MD');
ip.addRequired('refs');
ip.addParameter('name','projectFrameOfRef');
ip.addParameter('process',[])
ip.parse(MD,ref,varargin{:});
p=ip.Results;


%% Describe the Dynamical ROI ( A pyramid that descibe the maximum distance for randomization)
origTracks=TracksHandle();
nFrame=length(ref.frame);
origTracks.x=ones(1,nFrame);
origTracks.y=ones(1,nFrame);
origTracks.z=ones(1,nFrame);
origTracks.startFrame=1;
origTracks.endFrame=nFrame;
maxTracks=TracksHandle();
maxTracks.x=MD.imSize_(2)*ones(1,nFrame);
maxTracks.y=MD.imSize_(1)*ones(1,nFrame);
maxTracks.z=ones(1,nFrame*MD.zSize_*MD.pixelSizeZ_/MD.pixelSize_);
maxTracks.startFrame=1;
maxTracks.endFrame=nFrame;

dynFullPolygon=[ref.applyBase(origTracks,[]) ref.applyBase(maxTracks,[])];

project1D(  MD,[],'FoF',ref,'dynPoligonREF',dynFullPolygon, ...
    'name','fullSpindleRef','channelRender','grayRed', ...
    'processSingleProj',p.process,'intMinPrctil',[1 85],'intMaxPrctil',[100 100]);
