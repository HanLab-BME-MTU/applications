function [guiData] = ptRetrieveGUIData (handles)
% ptRetrieveGUIData loads data from the Polytrack_PP GUI 
%
% SYNOPSIS       [guiData, result] = ptRetrieveGUIData (handles)
%
% INPUT          handles : a pointer to the GUI data 
%
% OUTPUT         guiData : struct containing gui data related to jobs
%
% DEPENDENCIES   ptRetrieveGUIData uses {nothing}
%                                  
%                ptRetrieveGUIData is used by { PolyTrack_PP }
%
% Revision History
% Name                  Date            Comment
% --------------------- --------        --------------------------------------------------------
% Andre Kerstens        Sep 04          Initial Release
% Andre Kerstens        Mar 05          Added drug application timepoint

% Get the latest values from the GUI
guiData.plotfirstimg = str2num(get(handles.GUI_ad_firstimage_ed,'string'));
guiData.plotlastimg = str2num(get(handles.GUI_ad_lastimage_ed,'string'));
guiData.moviefirstimg = str2num(get(handles.GUI_fm_movieimgone_ed,'string'));
guiData.movielastimg = str2num(get(handles.GUI_fm_movieimgend_ed,'string'));
guiData.multframevelocity = str2num(get(handles.multFrameVelocity,'string'));
guiData.binsize = str2num(get(handles.GUI_ad_binsize_ed, 'String'));
guiData.dragtail = str2num(get(handles.GUI_fm_tracksince_ed, 'String'));
guiData.dragtailfile = get(handles.GUI_fm_filename_ed, 'String');
guiData.nrtrajectories = str2num(get(handles.nr_traj_ed, 'String'));
guiData.maxneighbourdist = str2num(get(handles.neighbour_dist_ed, 'String'));
guiData.maxcellcelldist = str2num(get(handles.GUI_maxcellcelldist_ed, 'String'));
guiData.drugtimepoint = str2num(get(handles.GUI_drugtimepoint_ed, 'String'));

% Get the path where to save new data
guiData.savedatapath = get(handles.GUI_fm_saveallpath_ed,'string');

% Get value of window size (for averaging)
guiData.windowsize = str2num(get(handles.GUI_windowsize_ed,'string'));

% Get value of confidence interval (for ripley)
guiData.ripleyconfint = str2num(get(handles.GUI_ripconfint_ed,'string'));

% Values for auto post processing
guiData.plusframes = str2num(get(handles.GUI_app_plusframes_ed, 'String'));
guiData.minusframes = str2num(get(handles.GUI_app_minusframes_ed, 'String'));
guiData.relinkdistance = str2num(get(handles.GUI_app_relinkdist_ed, 'String'));
guiData.mintrackdistance = str2num(get(handles.GUI_app_minimaltrack_ed, 'String'));
