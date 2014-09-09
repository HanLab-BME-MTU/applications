function [ filoInfo] = markFiloWithPersistentVeil(output,filoInfo,iFrame)
%UNTITLED9 Summary of this function goes here
%   Detailed explanation goes here
% for each filodopodia with a window figure out if that window is
% persistent (just add another field mark with a 1?) 
% the only thing is the output is currently for the whole movie (blockout)
% gives the frames of the window. 

% the easiest thing to do check if filo has a window, then check if that
% window is high persistence 
% get all filo indices (non-logical) that have an associated window. 
% branch structures are included in the filoopdia structure and will not
% have a window (ie only filoInfo.type <=1 will have a window)
idxBody = arrayfun(@(x) filoInfo(x).type<=1,1:length(filoInfo));
% have a small bug where some filo are not connected completely back to the
% base and hence are not assigned a window Idx filter and record here. 
idxWind = arrayfun(@(x) ~isempty(filoInfo(x).windowIdx),1:length(filoInfo)); 
idxNotWind = find(~idxWind);
% not sure how to assign all the data to structure 
for iNot = 1:length(idxNotWind)
    filoInfo(idxNotWind(iNot)).retractionpersVeil = 0; 
    filoInfo(idxNotWind(iNot)).protrusionpersVeil = 0; 
    
end 



if sum(idxWind) ~=0; 
    % write an error report.. 
    idxBad = find(idxBody & ~idxWind);
    arrayfun(@(x) display(['Filo ' num2str(idxBad(x)) 'in Frame' num2str(iFrame) 'is a body Filo but has no' ...
    'window']),1:length(idxBad)) ;
    
end 
idxWindFilo = find(idxBody & idxWind); % if is a veil/stem attached filo 
% should have a window- however sometimes I have a small bug that I need to
% rework. 20140909
% remember to set all the indices to zeros for these as well. 
persType{1} = 'retraction'; 
persType{2} = 'protrusion'; 
for i = 1:2
for ifilo = 1:length(idxWindFilo)
    % get the window associated with filo
    idxC = idxWindFilo(ifilo);
    fieldC = [persType{i} 'Analysis'];
    % get the time frames associated with that window that were undergoing 
    % persistent protrusion
    frameNumPers= output.(fieldC).windows(filoInfo(idxC).windowIdx).blockOut; 
    
    % extract from cells (cell number tells block idx - we don't care at
    % this point) 
    frameNumPers = vertcat(frameNumPers{:});
    currentFrame= sum(frameNumPers==iFrame) ; % think if there is a more effecient way to transfer over this 
    % info- the main problem is that the time frame in my structure is
    % embedded upstream in the analInfo not directly in the filoInfo. 
     
      fieldC2 = [persType{i} 'persVeil']; % note just thought it might be better 
        % to just associate with each filo the block stats? then can filter
        % here?.. this might be more intuitive though
    
    if currentFrame == 1
      
        filoInfo(idxC).(fieldC2) = 1; % the analInfo will encode the time frame. 
    else 
        filoInfo(idxC).(fieldC2) = 0 ; % no persistent veil at this time point
    end % think about if this is the best way to have this info.. as what is a 
    % persistent veil is going to change based on definition- need to have
    % this information more effectively marked etc to make user friendly. 

    
   


end % for ifilo
        
end  % for retract vs protrusion 
        
    end 
    
    



