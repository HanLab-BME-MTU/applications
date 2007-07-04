function makiCountSpots
%MAKICOUNTSPOTS counts spots for each frame of several movies
%
% SYNOPSIS: makiCountSpots
%
% INPUT 
%
% OUTPUT 
%
% REMARKS
%
% created with MATLAB ver.: 7.4.0.287 (R2007a) on Windows_NT
%
% created by: jdorn
% DATE: 03-Jul-2007
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% load data files
done = false;
data = struct('dataStruct',[]);
ct = 1;

while ~done
    data(ct).dataStruct = makiLoadDataFile;
    ans = questdlg('more data','','yes','no','cancel','yes');
    switch ans
        case 'yes'
            ct = ct + 1;
        case 'no'
            done = true;
        otherwise
            return
    end
end

% for every frame: count spots
nData = length(data);
for iData = 1:nData;
    % get nSpots
    nTimepoints = length(data(iData).dataStruct.initCoord);
    data(iData).nSpots = zeros(nTimepoints,1);
    for t = 1:nTimepoints
        data(iData).nSpots(t) = length(data(iData).dataStruct.initCoord{t});
    end
end

figure,
plot(data(1).nSpots,'xr'),
hold on,
plot(data(2).nSpots,'+g'),
plot(data(3).nSpots,'ob')



