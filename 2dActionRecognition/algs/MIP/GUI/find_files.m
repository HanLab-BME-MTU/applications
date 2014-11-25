function find_files(handles)
%Find movies file in start dir and present on GUI

    %Get inputs from GUI:
    params_file = getappdata(handles.GUIData,'params_file');
    run(params_file);
    start_dir = get(handles.start_dir,'String');
    
    %Find all movies (with suffix params.GUI.mov_suffix{}) recursively:
    files_list = search_files_r(start_dir,params.GUI.mov_suffix);
    %Convert to convinient DB representation:
    database_info.location = start_dir;
    short_names = cellfun(@(x)x(length(start_dir)+1:end),files_list,'UniformOutput',false);%remove long path
    rel_dir_names = cellfun(@(x)fileparts(x),short_names,'UniformOutput',false);%get relative dirs to start
    file_names = cellfun(@(x,y)x(length(y)+1:end),short_names,rel_dir_names,'UniformOutput',false);%only filenames
    %Create db_data:{filenames, rel_dirs}
    database_info.db_data = [file_names,rel_dir_names];
    
    %add results to GUI:
    setappdata(handles.GUIData,'database_info',database_info);
    set(handles.files_list,'String',short_names);
    set(handles.Status,'String',[num2str(length(short_names)),' files found.']);
end