function [data] = makeMockPattern(data, circleradius, circledist);
% make mock pattern and fill into the field .mockpattern
% 
% SYNOPSIS [data] = makeMockPattern(data, circleradius, circledist);
%
% INPUT     data:   experiment structure
%           circleradius:   radius (in pixel) of the pattern circles
%           circledist:     distance (in pixel) between circle centers
%
% OUTPUT    creates new field in the data structure called
%                .mockpattern  
%           which contains the mask image
%           
%
% Dinah Loerke, modified Sept 09, 2008

if ~isfield(data,'imagesize')
    data = determineImagesize(data);
end

% (relative) index positions of the pixels that have the specified distance
% from a pixel of position n in an image of this size
[miniImX, miniImY] = ndgrid(-circleradius:circleradius,-circleradius:circleradius);
miniDist = sqrt(miniImX.^2 + miniImY.^2);
%inside pixels
[XinPix, YinPix] = find(miniDist<=circleradius);
% subtract ceil(rad)+1 to express pixel positions as relative to a given
% position
XinPix = XinPix-(circleradius+1);
YinPix = YinPix-(circleradius+1);
% corresponding distances of all these points
DinPix = sqrt(XinPix.^2 + YinPix.^2);

% loop over all movies
for i=1:length(data)
    
    sx  = data(i).imagesize(1);
    sy  = data(i).imagesize(2);
    
    % initialize mask
    mockmask = zeros(sy,sx);
    
    % initialization positions
    init_posx = ceil(circledist*rand(1));
    init_posy = ceil(circledist*rand(1));
    
    % make grid of circle centers and shift by initialization positions
    [xgrid,ygrid] = ndgrid([0:circledist:sx], [0:circledist:sy]); 
    xgrid = xgrid + init_posx;
    ygrid = ygrid + init_posy;
    
    % loop over all circle centers
    for k=1:length(xgrid(:))
        % inside pixel values
        cinPix_x = XinPix + xgrid(k);
        cinPix_y = YinPix + ygrid(k);
        % determine positions that are also inside the image of the
        % specified size
        usepos = find( (cinPix_x>0) & (cinPix_y>0) & (cinPix_x<=sx) & (cinPix_y<=sy) );
        % set designated positions to 1
        for u=1:length(usepos)
            mockmask(cinPix_y(usepos(u)),cinPix_x(usepos(u))) = 1;
        end
    end
        
    data(i).mockmask = mockmask';
    
    
    
end

end % of function


%%========================================================================

function [m2]=DistanceMatrix(c1,c2);
%this subfunction makes a neighbour-distance matrix for input matrix c1
%(n1 x 2 points) and c2
%output: m2 (n1 x n1) matrix containing the distances of each point in c1 
%from each point in c2

[np1,sd1]=size(c1);
[np2,sd2]=size(c2);

m2=zeros(np1,np2);

for k = 1:np1
    for n = 1:np2
        d = sqrt((c1(k,1)-c2(n,1))^2+(c1(k,2)-c2(n,2))^2);
        m2(k,n)=d;
    end
end


end

        
    
