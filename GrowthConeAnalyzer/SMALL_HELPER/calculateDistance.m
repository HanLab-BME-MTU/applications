
function [dist,pixIdxTrunc,deltC] = calculateDistance(pixIdx,size,mkPlot,varargin)
% small function to calcuale the distance along a line of pixels (have to
% do this quite frequently so it was easier to make a small helper function ) 
% INPUT 
% pixIdx: 
% size :of image
% mkPlot : make a sanity check 
% getDistFromEnd : new param for the body estimation so that can get the
% medial axis transform distances for a body shape metric at a set length
% from the tip of the neurite 
% distCutOff: if getDistFromEnd is turned on requires the user input this
% cut-off 
% the order it spits out the coordinates is not correct (check this think
% it is fixed (20140812)

% OUTPUT 
% pixIdxTrunc: the input pixels cut-off from a certain length from position 1 
% this can then be used to collect the medial axis transform lengths 


%% 
ip = inputParser;
ip.addRequired('pixIdx',@(x) isvector(x) || isempty(x) );
ip.addRequired('size',@isvector);
ip.addRequired('mkPlot',@isscalar); % make optional later just keep this here now to keep consistency 
ip.addParamValue('distCutOff',[],@(x) isscalar(x) || isempty(x));
ip.addParamValue('pixelSize',0.216, @isscalar); % just make this default for me now 

ip.parse(pixIdx,size,mkPlot,varargin{:});



distCutOff = ip.Results.distCutOff; 
pixelSize = ip.Results.pixelSize; 




      [y,x]  = ind2sub(size,pixIdx);
     
     
     
       deltX = diff(x); 
       deltY = diff(y); 
    
    
    delt =  arrayfun(@(i) sqrt(deltX(i)^2+deltY(i)^2),1:length(deltX));
    deltC = cumsum(delt); 
    deltCMic = deltC.*pixelSize; 
    if ~isempty(distCutOff)
        
    
    pixIdxTrunc= pixIdx(deltCMic<distCutOff); 
    else 
        pixIdxTrunc = [];
    
    end 
    if ~isempty(delt) % ie not a singleton 
        
    dist = sum(delt); 
    else 
        dist = 0.001 ; % just set singleton distances to 0 for now : note the adjMat will not include any 
        % values of 0 so have to set to 0.001 
    end 
    dist = dist.*pixelSize; 
    if mkPlot ==1 
        
     plot(x,y,'color','r');
     scatter(x,y,10,'r','filled');
     text(x(1),y(1),num2str(dist),'color','g'); 
    end 
end 






