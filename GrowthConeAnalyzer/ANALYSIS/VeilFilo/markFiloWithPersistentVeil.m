function [ filoInfo] = markFiloWithPersistentVeil(output,filoInfo,iFrame)
%adds logical fields .protusionPersVeil and .retractionPersVeil assigning 
% each filo as either associated with persistent veil (defined in previous
% step) 1 or not 0. 
% 
% INPUT: output : is the output from getWindowsWithHighPersistence already defined 
%     what is a significant veil protrusion/retraction and segments out
%     blocks this function simply filters the blocks (not in the future maybe can 
%     just translate this info directly to the filoInfo and window? we'll
%     see.. 
%   filoInfo: Nx1 Structure array  associated with 1 frame with many fields describ
%            % where N = the number of filopodia structures measured in
%            that frame
%   iFrame: the frame associated with the filoInfo given - note this might be a 
%           bit of a shortcoming on how I currently set up my data
%           structures. filoInfo currently does not contain frame
%           information- that info is stored in the larger analInfo
%           structure. 
% for each filodopodia with a window figure out if that window is
% persistent (just add another field mark with a 1?) 
% the only thing is the output is currently for the whole movie (blockout)
% OUTPUT: filoInfo: with added fields

%%

% initiate field; 
for  i = 1:length(filoInfo) 
filoInfo(i).floatingFilo =0; % sorry for the loop just want logicals not empty so can easily cat without fucking with indexing. 
end 
% the easiest thing to do check if filo has a window, then check if that
% window is high persistence 
% get all filo indices (non-logical) that have an associated window. 
% branch structures are included in the filoopdia structure and will not
% have a window (ie only filoInfo.type <=1 will have a window)
% note idxBody currently for troubleshooting- a branch more towards the
% stem can also have a a window associated with it in the current
% implementation
 idxBody = arrayfun(@(x) filoInfo(x).type<=1,1:length(filoInfo)); 
% have a small bug where some filo are not connected completely back to the
% base and hence are not assigned a window Idx filter and record here. 

% get the indices of all filopodia associated with a window
idxWind = arrayfun(@(x) ~isempty(filoInfo(x).windowIdx),1:length(filoInfo)); 
idxNotWind = find(~idxWind); 

% set the 
for iNot = 1:length(idxNotWind)
    filoInfo(idxNotWind(iNot)).retractionpersVeil = 0; 
    filoInfo(idxNotWind(iNot)).protrusionpersVeil = 0; 
end 

% test for floaters
idxBad = find(idxBody & ~idxWind);

%ERROR REPORT FOR LATER : To REMOVE 20140911
if ~isempty(idxBad); 
    % write an error report.. 
    
    arrayfun(@(x) display(['Filo ' num2str(idxBad(x)) 'in Frame' num2str(iFrame) 'is a body Filo but has no' ...
    'window']),1:length(idxBad)) ;
for i = 1:length(idxBad) 
    filoInfo(idxBad(i)).floatingFilo = 1; 
end 
    %floatingFilo(countFloaters) = [iFrame,idxBad(x)]; 
    %countFloaters = countFloaters +1; 
end 

idxWindFilo = find(idxWind); % if is a veil/stem attached filo 
% again some branches may be close enough to the base to have an associated
% window- therefore do not want to filter by filo type
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
    
    



