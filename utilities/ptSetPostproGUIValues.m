function handles = ptSetPostproGUIValues (handles, nr)
% ptSetPostproGUIValues shows the current handles struct on the GUI
%
% SYNOPSIS       ptSetPostproGUIValues (handles)
%
% INPUT          handles : GUI handle struct containing the postpro info
%                nr : the number of the job to show
%
% OUTPUT         handles: updated handles struct
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
if handles.jobData(nr).imagesavailable == 1
    set (handles.GUI_pp_imagepath_ed, 'String', handles.jobData(nr).imagefilepath);
else
    set (handles.GUI_pp_imagepath_ed, 'String', '');
end
set (handles.pp_firstframe, 'String', handles.jobData(nr).firstimg);
set (handles.pp_lastframe, 'String', handles.jobData(nr).lastimg);
set (handles.pp_increment, 'String', handles.jobData(nr).increment);
set (handles.GUI_mmpixel_ed, 'String', handles.jobData(nr).mmpixel);
set (handles.GUI_frameinterval_ed, 'String', handles.jobData(nr).timeperframe);

set (handles.GUI_fm_saveallpath_ed, 'String', handles.guiData.savedatapath);
set (handles.GUI_ad_firstimage_ed, 'String', handles.guiData.plotfirstimg);
set (handles.GUI_ad_lastimage_ed, 'String', handles.guiData.plotlastimg);
set (handles.GUI_fm_movieimgone_ed, 'String', handles.guiData.moviefirstimg);
set (handles.GUI_fm_movieimgend_ed, 'String', handles.guiData.movielastimg);
set (handles.GUI_app_relinkdist_ed, 'String', handles.guiData.relinkdistance);
set (handles.GUI_app_minusframes_ed, 'String', handles.guiData.minusframes);
set (handles.GUI_app_plusframes_ed, 'String', handles.guiData.plusframes);
set (handles.GUI_app_minimaltrack_ed, 'String', handles.guiData.mintrackdistance);
set (handles.GUI_fm_tracksince_ed, 'String', handles.guiData.dragtaillength);
set (handles.GUI_fm_filename_ed, 'String', handles.guiData.dragtailfile);
set (handles.multFrameVelocity, 'String', handles.guiData.multframevelocity);
set (handles.GUI_ad_binsize_ed, 'String', handles.guiData.binsize);
set (handles.GUI_movietype_avi_rb, 'Value', 1);
set (handles.GUI_movietype_qt_rb, 'Value', 0);
set (handles.nr_traj_ed, 'String', handles.guiData.nrtrajectories);
set (handles.neighbour_dist_ed, 'String', handles.guiData.maxneighbourdist);
set (handles.GUI_windowsize_ed, 'String', handles.guiData.windowsize);
set (handles.GUI_maxcellcelldist_ed, 'String', handles.guiData.maxcellcelldist);
set (handles.GUI_ripconfint_ed, 'String', handles.guiData.ripleyconfint);
set (handles.pp_bad_frames,'String', handles.jobData(nr).nrbadframes);