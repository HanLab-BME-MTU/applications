% Example code for batch processing on protrusionAnalysis
% Last updated: August 26, 2008 by Shann-Ching Chen, LCCB

% flag for batch processing
handles.batch_processing = 1;

% mask directory
handles.directory_name = 'M:\hopkins\methodPaperwithHunter\cell4\masks';

% mask format
handles.FileType = '*.tif';

% result directory
handles.result_directory_name = 'M:\hopkins\methodPaperwithHunter\cell4\masks';

% time resolution, seconds per frame
handles.timevalue = 1;

% spatial resolution, nms per pixel
handles.resolutionvalue = 1;

% number of segments for activity map
handles.segvalue = 30;

% down sample to # of pixels if length for fragment is greater than #
handles.dl_rate = 10;


[OK, handles] = protrusionAnalysis(handles);
if ~OK
    errordlg(['Protrusion has not been done!'],'Error in Calculation','modal');
end

% results will be stored at 'M:\hopkins\methodPaperwithHunter\cell4\masks\analysis_dl10'
