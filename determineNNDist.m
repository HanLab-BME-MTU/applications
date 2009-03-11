function [data]=determineNNDist(data);
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


% number of entries in data
lens = length(data);

ordir = cd;

for i=1:lens
    
    
    fprintf(' movie %02d',i);   
    
    % number of frames for this exp
    lenf = data(i).movieLength;
    
    % current path
    currPath = data(i).source;
    cd(currPath);
    
    nndexist = 0;
    
    if exist('LifetimeInfo')==7
        cd('LifetimeInfo')
        if exist('NearNeiDistance.mat')==2
            nndexist = 1;
        end
        cd(currPath);
    end
    
    % if nearest neighbor info matrix doesn't already exist
    if nndexist == 0
        
        % check for TrackInfo
        if exist('TrackInfoMatrices')==7
            cd('TrackInfoMatrices')
            if exist('trackInfo.mat')==2
                % load trackInfo
                loadfile = load('trackInfo.mat');
                if isfield(loadfile,'trackInfo')
                    trackInfo = loadfile.trackInfo;
                elseif isfield(loadfile,'trackInfoMat')
                    trackInfo = loadfile.trackInfoMat;
                else
                    trackInfo = [];
                end        
            end
            cd(currPath);
        end   
            
        if ~isempty(trackInfo)
            % calculate nearest neighbor distance
            [NearNeiDistance] = calcNNDist(trackInfo);
            % make a directory where to save the data, unless it already exists
            if exist('LifetimeInfo')==7
                cd('LifetimeInfo');
            else
                mkdir('LifetimeInfo');
                cd('LifetimeInfo');
            end

            % save data
             save('NearNeiDistance.mat','NearNeiDistance');

        end % of if trackInfo isn't empty
    
    end % if lftInfo doesn't already exist
    
    cd(ordir);
    
    fprintf('\b\b\b\b\b\b\b\b\b');
    
end % of for i

fprintf('\n');

end % of function



%% ========================================================================
%
%                                SUBFUNCTIONS
%
% =========================================================================


function [nndistmat] = calcNNDist(trackinfo);
% for each point in each frame in tracking matrix, calculate nearest 
% neighbor distance
% INPUT: trackinfo matrix
% OUTPUT: nearest neighbor matrix

[stx,sty] = size(trackinfo);
tsiz = round(sty/8);
nearNeiMat = nan*zeros(stx,tsiz);

if (issparse(trackinfo))
    use_trackinfo = full(trackinfo);
else 
    use_trackinfo = trackinfo;
end

% loop over time points (columns)

for t=1:tsiz
    
    fprintf(' frame %04d',t);
    
    xi = 8*(t-1)+1;
    yi = 8*(t-1)+2;
    % find defined x coordinates
    xcoords = use_trackinfo(:,xi);
    xcoords_valpos = find(xcoords>0);
    uxcoords = xcoords(xcoords_valpos);
    % find defined y coordinates
    ycoords = use_trackinfo(:,yi);
    ycoords_valpos = find(ycoords>0);
    uycoords = ycoords(ycoords_valpos);
        
    % calculate distance matrix
    xyvec = [uxcoords uycoords];
    %[dm,vectorMatrix] = distMat(xyvec);
    [distMatrix]=distanceMatrixFunc(xyvec);
    
    distMatrix(distMatrix==0) = nan;
    distMatMin = nanmin(distMatrix,[],2);
    
    nndistmat(xcoords_valpos,t) = distMatMin;
    
   
    fprintf('\b\b\b\b\b\b\b\b\b\b\b');
    
end % of for loop over time points



end % of function



%==========================================================================


function[m2]=distanceMatrixFunc(c1,c2)
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
    
[ncx1,ncy1]=size(c1);
[ncx2,ncy2]=size(c2);
m2=zeros(ncx1,ncx2);

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