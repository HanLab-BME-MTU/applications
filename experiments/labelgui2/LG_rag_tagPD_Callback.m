function LG_rag_tagPD_Callback(hObject,eventdata,reAssignHandles)
%LG_rag_tagPD_Callback is the callback for the tag-pulldown-menus in reAssignGUI

% check the index of the current handle. If there is PD with the same
% index, switch its value with the previous value of the current handle

% find handles
pdHandles = reAssignHandles.pdHandles;

% index of current PD
myIdx = find(pdHandles == hObject);

% get all values
pdValues = LG_rag_getPdValues(pdHandles);

% find current value and mask it
myValue = pdValues(myIdx);

% don't swap if myValue is 'None'
if myValue ~= 1
    pdValues(myIdx) = -1; % pd-values can only be positive

    % find other PD with same value
    sameValueIdx = find(pdValues == myValue);

    % get old values
    oldPdValues = reAssignHandles.oldPdValues;

    % change selection to previous
    if ~isempty(sameValueIdx)
        set(pdHandles(sameValueIdx),'Value',oldPdValues(myIdx));
    end
end


% store current values as old values
reAssignHandles.oldPdValues = LG_rag_getPdValues(pdHandles);
guidata(hObject,reAssignHandles);
