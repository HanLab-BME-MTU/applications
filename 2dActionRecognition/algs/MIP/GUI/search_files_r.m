function files_list = search_files_r(start_path,suffix)
%find all files with suffix in current directory and add recursively next directories

    %stop when no inner dirs:
    cur_dir = dir(start_path);
    inner_map = ~(strcmp({cur_dir(:).name},'.') | strcmp({cur_dir(:).name},'..'));
    cur_dir = cur_dir(inner_map);
    
    %add my files, then inner dirs:
    files_list = {};
    for s=1:length(suffix)
        cur_suf_files = dir(fullfile(start_path,['*',suffix{s}]));
        files_list = [files_list;...
                      cellfun(@(x)fullfile(start_path,x),{cur_suf_files.name},'UniformOutput',false)'
                      ];
    end
    
    has_inner_dirs = any([cur_dir.isdir]);    
    if(~has_inner_dirs)
        %found leaf get out:
        return;
    end
    
    %go over inner dirs:
    inner_inds = find([cur_dir.isdir]);
    
    for d=1:length(inner_inds)
        new_list = search_files_r(fullfile(start_path,cur_dir(inner_inds(d)).name),suffix);
        files_list = [files_list;new_list];
    end
    
end