
function [coords] = getEndpoints(pixIdx,size,checkSingletons)
% checkSingletons input added for maintaining functionality with
% skel2Graph4OutGrowthMetric - here need to just get coordinates of
% singletons as well even though not technically an end point. 
%
if nargin<3 
    checkSingletons = 0; 
end 
     
singleFlag = 0 ;

if checkSingletons == 1
    if length(pixIdx) == 1
        singleFlag = 1;
        
    end
end

if singleFlag == 1
    
    % it is a singtlton
    [ye,xe]  = ind2sub(size,pixIdx);
else
    
    %
    endpoints = zeros(size);
    endpoints(pixIdx)=1;
    sumKernel = [1 1 1];
    % find endpoints of the floating candidates to attach (note in the
    % future might want to prune so that the closest end to the
    % body is the only one to be re-attatched: this will avoid double connections)
    endpoints = double((endpoints.* (conv2(sumKernel, sumKernel', padarrayXT(endpoints, [1 1]), 'valid')-1))==1);
    [ye,xe] = find(endpoints~=0);
end 
    coords(:,1) = xe;
    coords(:,2) = ye;
end



