[fileName,pathName] = uigetfile('Select the file contaning the nuclear marker channel');

[fileName,pathName] = uigetfile(bfGetFileExtensions, 'Select the file contaning the nuclear marker channel - histone-2B', handles.history.lastAnalyzedNucleusDir);
                                
dataFilePath = fullfile( pathName, fileName );

[imageData, metadata] = uiGetBioImageData(dataFilePath, 'channelDescriptionsAndNamesList', {{'Nuclear Marker Channel', 'channelIdNuclei'}});


% % pre-compute data needed for display
% handles = ComputeDisplayData( handles );
% 
% % Update handles structure
% guidata(hObject, handles);
% 
% % Run analysis
% RunAnalysis(hObject, handles);

% startMatlabPool();
