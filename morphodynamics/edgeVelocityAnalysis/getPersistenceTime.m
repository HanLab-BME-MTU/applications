function [Protrusion,Retraction,motionState] = getPersistenceTime(TS,deltaT,varargin)
%This function calculate the threshold to define a protrusion or retraction
%and uses it to estimate the time of protrusion and retraction for each
%event in TS
%
%IMPORTANT: This function has no time series pre-processing
%
%Usage:
%       [Protrusion,Retraction,Up,Dw] = getPersistenceTime(TS,deltaT,varargin)
%
% Input:
%       TS     - vector - Time series
%       deltaT - scalar - sampling rate (in seconds, every deltaT seconds)
%       Optional:
%               nStd - number of standard deviations used to calculate the
%               threshold
%
% Output:
%       Protrusion.PersTime  - vector - Protrusion time for all protrusive events
%       Protusion.BlockOut   - cell   - Time point where Protrusion happened
%       Protrusion.MaxVeloc  - max velocity within each protrusion block
%       Protrusion.Veloc     - mean velocity within each protrusion block 
%       Protrusion.MinVeloc  - min velocity within each protrusion block
%       Protrusion.MednVeloc - median velocity within each protrusion block
%       Protrusion.limit     - limit of the energy for the most fast component
% See also: findingProtRetrTime, getEdgeMotionPersistence
%
%Marco Vilela, 2012

%% Parsing the input ******************************************************
ip = inputParser;
ip.addRequired('TS',@isvector);
ip.addRequired('deltaT',@isscalar);
ip.addParamValue('nStd',1,@isscalar);
ip.addParamValue('plotYes',false,@islogical);

ip.parse(TS,deltaT,varargin{:});
nStd    = ip.Results.nStd;
plotYes = ip.Results.plotYes;
TS      = TS(:);

%**************************************************************************
%%
%Block of real numbers isolated by NaN
realBlock = findBlock(find(isfinite(TS)),5);% The constant 5 here defines the minimal number of points. This number is defined by the EMD algorithm. For more see  G. Rilling, P. Flandrin and P. Gonçalves"On Empirical Mode Decomposition and its algorithms",


gImf   = [];    
gNoise = [];
for iB   = 1:numel(realBlock)
    imfs = emd(TS(realBlock{iB}));
    gImf = [gImf imfs(1,:)];
    
    [~,noise] = testImf(imfs);
    gNoise    = [gNoise;noise];
end
%Mu       = mean(gImf);
sdtError = std(gNoise); 

%Defining the lower and upper noise confidence bands centered on 0
Protrusion.limit  =  sdtError*nStd; 
Retraction.limit  = -sdtError*nStd;
%*************************
 

TSprot  = NaN(size(TS));
TSretr  = NaN(size(TS));

TSprot(TS > Protrusion.limit) = TS(TS > Protrusion.limit);
TSretr(TS < Retraction.limit) = TS(TS < Retraction.limit);

ProtBlock = findBlock(find(isfinite(TSprot)),1 );
RetrBlock = findBlock(find(isfinite(TSretr)),1 );

Protrusion  = getStuff(Protrusion,ProtBlock,TS,deltaT);
Retraction  = getStuff(Retraction,RetrBlock,-TS,deltaT);

motionState = sum([TS > Protrusion.limit -(TS < Retraction.limit) ],2);

if plotYes
    
    figure
    plot(TS)
    hold on
    nProtB = numel(Protrusion.blockOut);
    nRetrB = numel(Retraction.blockOut);
    
    for i = 1:nProtB
        plot(Protrusion.blockOut{i},TS(Protrusion.blockOut{i}),'g','LineWidth',2)
    end
    for i = 1:nRetrB
        plot(Retraction.blockOut{i},TS(Retraction.blockOut{i}),'r','LineWidth',2)
    end
    
end

end%End of main function

function cellData =  getStuff(cellData,block,TS,deltaT)

if ~isempty(block)
    
    [aux] = findingProtRetrTime(block,TS,deltaT);
    cellData = mergestruct(cellData,aux);
else
    
    cellData.persTime  = NaN;
    cellData.blockOut  = {[]};
    cellData.maxVeloc  = NaN;
    cellData.Veloc     = NaN;
    cellData.minVeloc  = NaN;
    cellData.mednVeloc = NaN;
    cellData.maxTime   = NaN;
end

end