function newSpeckleArray=fsmConvertSpeckleArray(speckleArray)
% fsmConvertSpeckleArray converts the speckleArray structure to the new format
%
% The old speckleArray was an 'array of structs'; the new one is a 'struct of arrays', 
% a format which guarantees a huge reduction in memory usage.
%
% SYNOPSIS      newSpeckleArray=fsmConvertSpeckleArray(speckleArray)
%
% INPUT         speckleArray    : speckle structure (see fsmBuildSaveSpeckleArray for more info)
%
% OUTPUT        newSpeckleArray : speckleArray with the new format
%
% DEPENDENCES   fsmConvertSpeckleArray uses { }
%               fsmMain is used by { fsmCenter }
%
% Aaron Ponti, May 7th, 2004

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Check input
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin~=1
    error('One input parameter expected.');
end

% Check for fields
fields={'timepoint','spPos','bgPos1','bgPos2','bgPos3','intensity', ...
        'background','deltaI','deltaICrit','sigmaSp','sigmaBg', ...
        'status','speckleType','score','activity','lmEvent'};
for i=1:length(fields)
    if isfield(speckleArray,char(fields(i)))==0
        error('The input speckleArray is not valid.');
    end
end

% Check that the speckleArray is really in the old format
if length(speckleArray)==1 & length(speckleArray(1).timepoint)>1
    error('This speckleArray is already in the new format.');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Initialize empty speckle array
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

lengthSpeckleArray=length(speckleArray);
newSpeckleArray=struct(...
    'timepoint',   uint16(zeros(lengthSpeckleArray,1)),...  % Unsigned integer 16: 1 - 65535 (2 bytes)
    'spPos',       uint16(zeros(lengthSpeckleArray,2)),...  % Unsigned integer 16 
    'bgPos1',      uint16(zeros(lengthSpeckleArray,2)),...  % Unsigned integer 16
    'bgPos2',      uint16(zeros(lengthSpeckleArray,2)),...  % Unsigned integer 16
    'bgPos3',      uint16(zeros(lengthSpeckleArray,2)),...  % Unsigned integer 16
    'intensity',   zeros(lengthSpeckleArray,1),...          % Double (8 bytes)
    'background',  zeros(lengthSpeckleArray,1),...          % Double
    'deltaI',      zeros(lengthSpeckleArray,1),...          % Double
    'deltaICrit',  zeros(lengthSpeckleArray,1),...          % Double
    'sigmaSp',     zeros(lengthSpeckleArray,1),...          % Double 
    'sigmaBg',     zeros(lengthSpeckleArray,1),...          % Double
    'status',      char(zeros(lengthSpeckleArray,1)),...    % Char (1 byte)
    'speckleType', uint8(zeros(lengthSpeckleArray,1)),...   % Unsigned integer 8: 1 - 255 (1 byte)
    'score',       zeros(lengthSpeckleArray,1),...          % Double 
    'activity',    int8(zeros(lengthSpeckleArray,1)),...    % Integer 8: -128 - 127 (1 byte)
    'lmEvent',     logical(zeros(lengthSpeckleArray,1)));   % Logical: 0 | 1 (1 byte)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Import fields
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

newSpeckleArray.timepoint(1:lengthSpeckleArray)   = [speckleArray.timepoint];
newSpeckleArray.spPos(1:lengthSpeckleArray,1:2)   = reshape([speckleArray.spPos],2,lengthSpeckleArray)';
newSpeckleArray.bgPos1(1:lengthSpeckleArray,1:2)  = reshape([speckleArray.bgPos1],2,lengthSpeckleArray)';
newSpeckleArray.bgPos2(1:lengthSpeckleArray,1:2)  = reshape([speckleArray.bgPos2],2,lengthSpeckleArray)';
newSpeckleArray.bgPos3(1:lengthSpeckleArray,1:2)  = reshape([speckleArray.bgPos3],2,lengthSpeckleArray)';
newSpeckleArray.intensity(1:lengthSpeckleArray)   = [speckleArray.intensity];
newSpeckleArray.background(1:lengthSpeckleArray)  = [speckleArray.background];
newSpeckleArray.deltaI(1:lengthSpeckleArray)      = [speckleArray.deltaI];
newSpeckleArray.deltaICrit(1:lengthSpeckleArray)  = [speckleArray.deltaICrit];
newSpeckleArray.sigmaSp(1:lengthSpeckleArray)     = [speckleArray.sigmaSp];
newSpeckleArray.sigmaBg(1:lengthSpeckleArray)     = [speckleArray.sigmaBg];
newSpeckleArray.status(1:lengthSpeckleArray)      = [speckleArray.status];
newSpeckleArray.speckleType(1:lengthSpeckleArray) = [speckleArray.speckleType];
newSpeckleArray.score(1:lengthSpeckleArray)       = [speckleArray.score];
newSpeckleArray.activity(1:lengthSpeckleArray)    = [speckleArray.activity];
newSpeckleArray.lmEvent(1:lengthSpeckleArray)     = [speckleArray.lmEvent];

