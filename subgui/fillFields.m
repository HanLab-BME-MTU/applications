function fillFields(handles,activeJob)


set(handles.GUI_st_path_imagedirectory_ed,'String',activeJob.imagedirectory);
set(handles.GUI_st_path_imagename_ed,'String',activeJob.imagename);
set(handles.GUI_st_path_firstimage_ed,'String',num2str(activeJob.firstimage));
set(handles.GUI_st_path_lastimage_ed,'String',num2str(activeJob.lastimage));
set(handles.GUI_st_path_increment_ed,'String',num2str(activeJob.increment));
set(handles.GUI_st_path_savedirectory_ed,'String',activeJob.savedirectory);
set(handles.GUI_st_iq_fi_background_ed,'String',num2str(activeJob.fi_background));
set(handles.GUI_st_iq_fi_nucleus_ed,'String',num2str(activeJob.fi_nucleus));
set(handles.GUI_st_iq_la_background_ed,'String',num2str(activeJob.la_background));
set(handles.GUI_st_iq_la_nucleus_ed,'String',num2str(activeJob.la_nucleus));
set(handles.GUI_st_bp_maxsearch_ed,'String',num2str(activeJob.maxsearch));
set(handles.GUI_st_bp_minsize_ed,'String',num2str(activeJob.minsize));
set(handles.GUI_st_bp_maxsize_ed,'String',num2str(activeJob.maxsize));
set(handles.GUI_st_bp_minsdist_ed,'String',num2str(activeJob.minsdist));
set(handles.GUI_st_bp_halolevel_ed,'String',num2str(activeJob.halolevel));
set(handles.GUI_st_eo_minedge_ed,'String',num2str(activeJob.minedge));
set(handles.GUI_st_eo_sizetemplate_ed,'String',num2str(activeJob.sizetemplate));
set(handles.GUI_st_eo_mintrackcorrqual_ed,'String',num2str(activeJob.mintrackcorrqual));
set(handles.GUI_st_eo_boxsize_ed,'String',num2str(activeJob.boxsize));

set(handles.GUI_st_eo_mincorrqualtempl_pm,'String',num2str(activeJob.mincorrqualtempl));
set(handles.GUI_st_eo_noiseparameter_pm,'String',num2str(activeJob.noiseparameter));
set(handles.GUI_st_eo_leveladjust_pm,'String',num2str(activeJob.leveladjust));



if ~isempty(activeJob.mmpixel)
set(handles.GUI_st_bp_mmpixel_pm,'String',num2str(activeJob.mmpixel));
end

if ~isempty(activeJob.timestepslide)
set(handles.GUI_st_eo_timestepslide_pm,'String',num2str(activeJob.timestepslide));
end



