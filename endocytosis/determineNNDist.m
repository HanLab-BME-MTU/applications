function [data] = determineNNDist(data)
% determine nearest neighbor distance for all entries in the data structure
%
% SYNOPSIS  [data]=determineNNDist(data);
%
% INPUT     data    :    experiment structure, which has to contain a field
%                       .source, which is the path to the data location; at
%                       this location, the function reads the trackInfo
%                       from a folder called TrackInfoMatrices
% OUTPUT    data           
% REMARKS 
%
% Dinah Loerke, last modified Feb 2008
% Francois Aguet, Jan 2010

for i = 1:length(data)
    fprintf(' movie %02d',i);
    
    % if nearest neighbor info matrix doesn't already exist
    if ~(exist([data(i).source 'LifetimeInfo' filesep 'NearNeiDistance.mat'], 'file')==2)
        
        trackInfoPath = [data(i).source 'TrackInfoMatrices' filesep 'trackInfo.mat'];
        if (exist(trackInfoPath, 'file')==2)
            loadfile = load(trackInfoPath);
            if isfield(loadfile,'trackInfo')
                NearNeiDistance = calcNNDist(loadfile.trackInfo);
            elseif isfield(loadfile,'trackInfoMat')
                NearNeiDistance = calcNNDist(loadfile.trackInfoMat);
            else
                error('No suitable trackInfo.mat file found.');
            end
            if ~(exist([data(i).source 'LifetimeInfo'], 'dir')==7)
                mkdir([data(i).source 'LifetimeInfo']);
            end
            save([data(i).source 'LifetimeInfo' filesep 'NearNeiDistance.mat'], 'NearNeiDistance')
        end
    end
end
end % of function


function [nndistmat] = calcNNDist(trackinfo)
% for each point in each frame in tracking matrix, calculate nearest 
% neighbor distance
% INPUT: trackinfo matrix
% OUTPUT: nearest neighbor matrix

if (issparse(trackinfo))
    use_trackinfo = full(trackinfo);
else 
    use_trackinfo = trackinfo;
end

% loop over time points (columns)

for t = 1:round(size(trackinfo, 2)/8)
    
    fprintf(' frame %04d',t);
    
    xi = 8*(t-1)+1;
    yi = 8*(t-1)+2;
    % find defined x coordinates
    xcoords = use_trackinfo(:,xi);
    xcoords_valpos = find(xcoords>0);
    uxcoords = xcoords(xcoords_valpos);
    % find defined y coordinates
    ycoords = use_trackinfo(:,yi);
    uycoords = ycoords(ycoords>0);
        
    % calculate distance matrix
    xyvec = [uxcoords uycoords];
    %[dm,vectorMatrix] = distMat(xyvec);
    [distMatrix]=distanceMatrixFunc(xyvec);
    
    distMatrix(distMatrix==0) = nan;
    distMatMin = nanmin(distMatrix,[],2);
    
    nndistmat(xcoords_valpos,t) = distMatMin;
    
   
    fprintf('\b\b\b\b\b\b\b\b\b\b\b');
    
end
end % of function



%==========================================================================


function[m2] = distanceMatrixFunc(c1, c2)
%this subfunction makes a neighbour-distance matrix for input matrix m1
%input: c1 (n1 x 2 points) and c2 (n2 x 2 points) matrices containing 
%the x,y coordinates of n1 or n2 points
%output: m2 (n1 x n2) matrix containing the distances of each point in c1 
%from each point in c2

vself = 0;
if nargin<2
    c2 = c1;
    vself = 1;
end
    
ncx1 = size(c1,1);
ncx2 = size(c2,1);
m2 = zeros(ncx1,ncx2);

if vself == 0
    for k=1:ncx1
        for n=1:ncx2
            d=sqrt((c1(k,1)-c2(n,1))^2+(c1(k,2)-c2(n,2))^2);
            m2(k,n)=d;
        end
    end
else
    for k=1:ncx1
        for n=k+1:ncx2
            d=sqrt((c1(k,1)-c2(n,1))^2+(c1(k,2)-c2(n,2))^2);
            m2(k,n)=d;
            m2(n,k)=d;
        end
    end
end
    
end % of subfunction