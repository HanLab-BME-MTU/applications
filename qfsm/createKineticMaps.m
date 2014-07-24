function [polyMap,depolyMap,kinMap2C]=createKineticMaps(kinScore,n,imSize,varargin)
% fsmKineticMaps creates maps of polymerization, depolymerization and net assembly rate
%
% SYNOPSIS      [polyMap,depolyMap,kinMap2C]=fsmKineticMaps(projDir,n,sigma)
%
% INPUT         kinScore        : string pointing to the project directory (where the fsmParam.mat file is 
%                                located). Pass projDir=[] to manually select fsmParam.mat via a dialog.
%               n              : number of frames for time integration
%               imgSize   : pixel size in the image domain (nm)
%               sigma (pixels) : sigma for the low-pass filtering of the maps
%                                (default - sigma=5)
%
% OUTPUT        polyMap        : 2D polymerization map integrated over n frames
%               depolyMap      : 2D depolymerization map integrated over n frames
%               netMap         : dual color 2D kinetic map integrated over n frames
%                                red channel for polymerization; green channel for depolymerization
%                                The netMap maps are always normalized with respect to the highest 
%                                score (positive or negative); polyMap and depolyMap are not stretched
%
% DEPENDENCES   fsmKineticMaps uses { }
%               fsmKineticMaps is used by {  }
%
% Aaron Ponti, September 2th, 2003
% Sebastien Besson, June 2011
% Adapted from fsmKineticMaps

ip = inputParser;
ip.addRequired('kinScore',@iscell);
ip.addRequired('n',@(x) isscalar(x) | mod(x,2)~=0);
ip.addRequired('imSize',@isvector);
ip.addOptional('sigma',5,@isscalar);
ip.addParamValue('waitbar',[],@ishandle)
ip.parse(kinScore,n,imSize,varargin{:});
sigma = ip.Results.sigma;

if ~isempty(ip.Results.waitbar)
    wtBar=ip.Results.waitbar;
    waitbar(0,wtBar,'Creating kinetic maps');
elseif feature('ShowFigureWindows')
    wtBar = waitbar(0,'Creating kinetic maps');
end

% Initialize the output
nMaps=numel(kinScore)-n+1;
polyMap=repmat({zeros(imSize)},1,nMaps);
depolyMap=repmat({zeros(imSize)},1,nMaps);
kinMap2C=repmat({zeros([imSize 3])},1,nMaps);


% Find all polymerization and depolymerization events
polyEvents = cellfun(@(x) x(x(:,4)>0,:),kinScore,'UniformOutput',false);
depolyEvents = cellfun(@(x) x(x(:,4)<0,:),kinScore,'UniformOutput',false);

for i=1:nMaps
    
    % Create polymerization map
    polyInd=cellfun(@(x)sub2ind(imSize,x(:,2),x(:,3)),polyEvents(i:i+n-1),...
        'UniformOutput',false);
    polyScore=cellfun(@(x) x(:,4),polyEvents(i:i+n-1),'UniformOutput',false);
    for j=1:n
        polyMap{i}(polyInd{j})=polyMap{i}(polyInd{j})+polyScore{j};
    end
    
    % Create depolymerization map
    depolyInd=cellfun(@(x)sub2ind(imSize,x(:,2),x(:,3)),depolyEvents(i:i+n-1),...
        'UniformOutput',false);
    depolyScore=cellfun(@(x) x(:,4),depolyEvents(i:i+n-1),'UniformOutput',false);
    for j=1:n
        depolyMap{i}(depolyInd{j})=depolyMap{i}(depolyInd{j})+depolyScore{j};
    end
    
    % Average maps
    polyMap{i}=polyMap{i}/n;
    depolyMap{i}=depolyMap{i}/n;
    
    % Low-pass filter
    if sigma~=0
        polyMap{i}=filterGauss2D(polyMap{i},sigma);
        depolyMap{i}=filterGauss2D(depolyMap{i},sigma);
    end
    
    % Create dual-channel image
    mx=max([polyMap{i}(:);abs(depolyMap{i}(:))]);
    kinMap2C{i}(:,:,1)=polyMap{i}/mx;
    kinMap2C{i}(:,:,2)=abs(depolyMap{i})/mx;
    
    % Use NaN in mpas for undetected events
    polyMap{i}(polyMap{i} == 0) = NaN;
    depolyMap{i}(depolyMap{i} == 0) = NaN;
    
    % Update waitbar if applicable
    if mod(i,round(nMaps/20))==1 && ishandle(wtBar),
        waitbar(i/nMaps,wtBar); 
    end
end

% Close waitbar if not-delegated
if isempty(ip.Results.waitbar) && ishandle(wtBar), 
    close(wtBar); 
end