function varargout = curveDirection3D(polyXYZ,method,projDist)
%POLYDIRECTION3D gets average / approxiamte direction of a 3d curve
%
%[dirVect,projLine] = polyDirection3D(polyXYZ,method,projDist) 
%More HELP coming sooon.....
%

if nargin < 1 || isempty(polyXYZ) || size(polyXYZ,2) ~= 3 || size(polyXYZ,1) <2
    error('The first input must be an Mx3 matrix, where M is >= 2!')
end

if nargin < 2 || isempty(method)
    method = 'Average';
end

if nargin < 3 || isempty(projDist)
    projDist = 10;
end

nPts = size(polyXYZ,1);

showPlots = false;

switch method
    
    
    case 'Average'

        %Get averaged direction vector
        dirVect = mean(diff(polyXYZ,1,1),1);
        
    case 'WeightedAverage'
        
        %Get average direction, weighting later points more highly.        
        
        %Create matrix for weighting points
        wtMat = repmat((1:nPts-1)' .^2,[1 3]);
        dirVect = sum(diff(polyXYZ,1,1) .* wtMat,1) ./ sum(wtMat,1);
        
    otherwise
        
        if isposint(method)
        
            %Get average of last few points
            nMean = min(method,nPts-1);
            dirVect = mean(diff(polyXYZ(end-nMean:end,:),1,1),1);
        else                    
            error(['"' method '" is not a recognized method!'])
        end
end


dirVect = dirVect ./ norm(dirVect);

if nargout > 0
    varargout{1} = dirVect;
end
if nargout > 1
    
    nProj = ceil(projDist);
    
    %Get projected points in this direction with unit spacing
    projLine = repmat(polyXYZ(end,:),[nProj+1 1])  + repmat(dirVect,[nProj+1 1]) .* repmat((0:nProj)',[1 3]);
    varargout{2} = projLine;
    
    
end

        
if showPlots
    
    %Get first-to-last displacement vector
    dispVect = diff(polyXYZ([1,end],:),1,1);

    vecOut = (dirVect ./ norm(dirVect)) .* norm(dispVect);
       
    figure;
    plot3(polyXYZ(:,1),polyXYZ(:,2),polyXYZ(:,3),'.-b')
    hold on    
    quiver3(polyXYZ(end,1),polyXYZ(end,2),polyXYZ(end,3),vecOut(1),vecOut(2),vecOut(3),0,'r')
    if nargout > 1
        plot3(projLine(:,1),projLine(:,2),projLine(:,3),'.-g')
    end
    axis equal
    axis vis3d
    
end