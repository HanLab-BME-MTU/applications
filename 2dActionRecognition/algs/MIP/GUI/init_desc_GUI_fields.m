function init_desc_GUI_fields(handles)
%Init GUI fields.
    %get params:
    params_file = getappdata(handles.GUIData,'params_file');
    run(params_file);
    
    %Set supported datasets acourding to the params file:
    cur_str = get(handles.datasets_select,'String');
    new_str = params.GUI.supported_datasets(:,1);%names only    
    val = get(handles.datasets_select,'Value');
    %set the new value:
    if(~isempty(cur_str) && iscell(cur_str))
        cur_selected = cur_str{val};
    else
        cur_selected = cur_str;
    end
    new_val = find(strcmp(new_str,cur_selected));
    if(isempty(new_val) || new_val < 1)
        new_val = 1;
    end
    set(handles.datasets_select,'String',new_str);
    set(handles.datasets_select,'Value',new_val);
    update_dataset_select(handles);


end