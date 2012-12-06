function [Protrusion,Retraction] = getPersistenceTime(TS,deltaT,varargin)
%This function calculate the threshold to define a protrusion or retraction
%and uses it to estimate the time of protrusion and retraction for each
%event in TS
%Usage:
%       [Protrusion,Retraction,Up,Dw] = getPersistenceTime(TS,deltaT,varargin)
%
% Input:
%       TS     - vector - Time series
%       deltaT - scalar - sampling rate (in seconds, every deltaT seconds)
%       Optional:
%               per - number of standard deviations used to calculate the
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
ip.addOptional('per',1,@isscalar);

ip.parse(TS,deltaT,varargin{:});
per  = ip.Results.per;
%**************************************************************************
%%

imf       = emd(TS);
Mu        = mean(imf(1,:));
nPoint    = length(TS);
[~,noise] = testImf(imf);
sdtError  = std(noise); 

%Defining the lower and upper noise confidence bands
Protrusion.limit  = Mu + sdtError*per; 
Retraction.limit  = Mu - sdtError*per;
%*************************
 

TSprot  = NaN(size(TS));
TSretr  = NaN(size(TS));

TSprot(TS > Protrusion.limit) = TS(TS > Protrusion.limit);
TSretr(TS < Retraction.limit) = TS(TS < Retraction.limit);

ProtBlock = findBlock( setdiff(1:nPoint,find(isnan(TSprot))),1 );
RetrBlock = findBlock( setdiff(1:nPoint,find(isnan(TSretr))),1 );

Protrusion = getStuff(Protrusion,ProtBlock,TS,deltaT);
Retraction = getStuff(Retraction,RetrBlock,-TS,deltaT);

end%End of main function

function cellData =  getStuff(cellData,block,TS,deltaT)

if ~isempty(block)
    
    [cellData.persTime,cellData.blockOut,cellData.maxVeloc,cellData.Veloc,cellData.minVeloc,cellData.mednVeloc] ...
            = findingProtRetrTime(block,TS,deltaT);
    
else
    
    cellData.persTime  = NaN;
    cellData.blockOut  = {[]};
    cellData.maxVeloc  = NaN;
    cellData.Veloc     = NaN;
    cellData.minVeloc  = NaN;
    cellData.mednVeloc = NaN;
    
end

end