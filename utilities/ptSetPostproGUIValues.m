function ptSetPostproGUIValues (handles)
% ptSetPostproGUIValues shows the current handles struct on the GUI
%
% SYNOPSIS       ptSetPostproGUIValues (handles)
%
% INPUT          handles : GUI handle struct containing the postpro info
%
% OUTPUT         None
%
% DEPENDENCIES   ptSetPostproGUIValues uses {nothing}
%                                  
%                ptSetPostproGUIValues is used by { PolyTrack_PP }
%
% Revision History
% Name                  Date            Comment
% --------------------- --------        --------------------------------------------------------
% Andre Kerstens        Aug 04          Initial Release

% Update fields on the GUI with the latest values
%set (handles.GUI_pp_jobpath_ed, 'String', jobValPath);
set (handles.GUI_pp_imagepath_ed, 'String', handles.postpro.imagepath);
set (handles.GUI_fm_saveallpath_ed, 'String', handles.saveallpath);
set (handles.GUI_ad_firstimage_ed, 'String', handles.postpro.plotfirstimg);
set (handles.GUI_ad_lastimage_ed, 'String', handles.postpro.plotlastimg);
set (handles.GUI_fm_movieimgone_ed, 'String', handles.postpro.moviefirstimg);
set (handles.GUI_fm_movieimgend_ed, 'String', handles.postpro.movielastimg);
set (handles.GUI_app_relinkdist_ed, 'String', handles.postpro.maxdistpostpro);
set (handles.GUI_app_minusframes_ed, 'String', handles.postpro.minusframes);
set (handles.GUI_app_plusframes_ed, 'String', handles.postpro.plusframes);
set (handles.GUI_app_minimaltrack_ed, 'String', handles.postpro.minimaltrack);
set (handles.GUI_fm_tracksince_ed, 'String', handles.postpro.dragtail);
set (handles.GUI_fm_filename_ed, 'String', handles.postpro.dragtailfile);
set (handles.multFrameVelocity, 'String', handles.postpro.multFrameVelocity);
set (handles.GUI_ad_binsize_ed, 'String', handles.postpro.binsize);
set (handles.pp_firstframe, 'String', handles.postpro.firstimg);
set (handles.pp_lastframe, 'String', handles.postpro.lastimg);
set (handles.pp_increment, 'String', handles.postpro.increment);
set (handles.GUI_mmpixel_ed, 'String', handles.postpro.mmpixel);
set (handles.GUI_frameinterval_ed, 'String', handles.postpro.timeperframe);
set (handles.GUI_movietype_avi_rb, 'Value', 1);
set (handles.GUI_movietype_qt_rb, 'Value', 0);
set (handles.nr_traj_ed, 'String', handles.postpro.nrtrajectories);
set (handles.neighbour_dist_ed, 'String', handles.postpro.neighbourdist);
set (handles.GUI_windowsize_ed, 'String', handles.postpro.windowsize);
