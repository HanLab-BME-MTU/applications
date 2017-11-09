function tracksOut = trackScalar(tracks,scalar)
%TRACKSCALAR Multiply track coordinates by a scalar
%   tracksOut = trackScalar(tracks,scalar)
%
%   Inputs:
%       tracks:     A track structure containing the field
%                   'tracksCoordAmpCG'
%
%       scalar:     A constant by which to multiply all values in
%                   tracks.tracksCoordAmpCG
%
%   Output:
%       tracksOut:  New track structure with scaled coordinates
%
%Kevin Nguyen, July 2016

% Multiply scalar to every compound track 
newCoord = arrayfun(@(x) x.tracksCoordAmpCG*scalar,tracks,'UniformOutput',false);
% newCoordStruct = cell2struct(newCoord','tracksCoordAmpCG',1);

tracksOut = tracks;
[tracksOut.tracksCoordAmpCG] = newCoord{:};

end