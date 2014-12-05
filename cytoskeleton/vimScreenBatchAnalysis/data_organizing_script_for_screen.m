% Data organizing script for screen data

% Liya Ding 2014.12

%Define directories

Data_Dir = '/project/cellbiology/gdanuser/vimentin/tonyzhang/live cell nikon data/well plate screening/screen_20141204';
%
Analysis_Dir = '/project/cellbiology/gdanuser/vimentin/tonyzhang/live cell nikon data/well plate screening/screen_20141204_analysis';

% Data_Dir = 'D:\Work\fromTony_201412';
%
% Analysis_Dir = 'D:\Work\fromTony_201412_analysis';

if(~exist(Analysis_Dir,'dir'))
    mkdir(Analysis_Dir);
end

%% First, script for changing the folder names with zero padding

% the well names are A2-Site_2 or B10-Site_3
% rename them to have common form so that A2-Site_2 will be A02-Site_2

data_folder_list =  dir([Data_Dir,filesep,'*Site*']);

data_folder_number = numel(data_folder_list);

correct_data_folder_list = cell(data_folder_number,1);

for iD = 1 : data_folder_number
    folder_name = data_folder_list(iD).name;
    correct_data_folder_list{iD} =  folder_name;
    
    % for column index and missing 0 paddings
    ind = find(folder_name=='-');
    
    if(ind(1)==3)
        folder_name = [folder_name(1) '0' folder_name(2:end)];
        correct_data_folder_list{iD} =  folder_name;
    end
end

%% Copy and Rename the images of three channels from different folders
% to three big folders, organized by each row.
% So all row A will be in one big folder and one movieData object.
% and there will be 16 rows, so 16 movieData object into a movieList for a
% plate.

% now corrected with 0 padding

for iD = 1 : data_folder_number
    folder_name = correct_data_folder_list{iD};
    
    Row_Name(iD,1) = folder_name(1);
    
    if(~exist([Analysis_Dir,filesep,Row_Name(iD,1),'_Row_Analysis'],'dir'))
        mkdir([Analysis_Dir,filesep,Row_Name(iD,1),'_Row_Analysis']);
    end
    
    if(~exist([Analysis_Dir,filesep,Row_Name(iD,1),'_Row_Analysis',filesep,'Vim_Images'],'dir'))
        mkdir([Analysis_Dir,filesep,Row_Name(iD,1),'_Row_Analysis',filesep,'Vim_Images']);
    end
    
    if(~exist([Analysis_Dir,filesep,Row_Name(iD,1),'_Row_Analysis',filesep,'MT_Images'],'dir'))
        mkdir([Analysis_Dir,filesep,Row_Name(iD,1),'_Row_Analysis',filesep,'MT_Images']);
        
    end
    
    if(~exist([Analysis_Dir,filesep,Row_Name(iD,1),'_Row_Analysis',filesep,'DAPI_Images'],'dir'))
        mkdir([Analysis_Dir,filesep,Row_Name(iD,1),'_Row_Analysis',filesep,'DAPI_Images']);
    end
    
    % already added 0 padding, so always 2:3
    Column_Name = folder_name(2:3);
    Column_Ind(iD,1) = str2num(Column_Name);
    
    % for site index
    Site_Name = folder_name(end);
    Site_Ind(iD,1) = str2num(Site_Name);
    
    % copy only if the three files are all present
    if(exist([Data_Dir, filesep,data_folder_list(iD).name, filesep,'img_000000000_FITC_000.tif'],'file'))
        if(exist([Data_Dir, filesep,data_folder_list(iD).name, filesep,'img_000000000_TRITC_000.tif'],'file'))
            if(exist([Data_Dir, filesep,data_folder_list(iD).name, filesep,'img_000000000_Dapi_000.tif'],'file'))
                
                copyfile([Data_Dir, filesep,data_folder_list(iD).name, filesep,'img_000000000_FITC_000.tif'], ...
                    [Analysis_Dir,filesep,Row_Name(iD,1),'_Row_Analysis',filesep,'Vim_Images',filesep,'vim_',Row_Name(iD,1),'_',Column_Name,'_s',Site_Name,'.tif']);
                
                copyfile([Data_Dir, filesep,data_folder_list(iD).name, filesep,'img_000000000_TRITC_000.tif'], ...
                    [Analysis_Dir,filesep,Row_Name(iD,1),'_Row_Analysis',filesep,'MT_Images',filesep,'mt_',Row_Name(iD,1),'_',Column_Name,'_s',Site_Name,'.tif']);
                
                copyfile([Data_Dir, filesep,data_folder_list(iD).name, filesep,'img_000000000_Dapi_000.tif'], ...
                    [Analysis_Dir,filesep,Row_Name(iD,1),'_Row_Analysis',filesep,'DAPI_Images',filesep,'dapi_',Row_Name(iD,1),'_',Column_Name,'_s',Site_Name,'.tif']);
            end
        end
    end
end


%% Build the movieData object and movieList for a  plate.

all_rows = unique(Row_Name);

Movie_number = numel(all_rows);

movie_cell = cell(1,Movie_number);

for iMovie = 1 : Movie_number
    
    %% Add the three channels
    channels_obj_all=[];
    
    channels_obj_one = Channel([Analysis_Dir,filesep,all_rows(iM),'_Row_Analysis',filesep,'Vim_Images']);
    channels_obj_one.getImageFileNames();
    channels_obj_one.sanityCheck();
    channels_obj_all = [channels_obj_all channels_obj_one];
    
    channels_obj_one = Channel([Analysis_Dir,filesep,all_rows(iM),'_Row_Analysis',filesep,'MT_Images']);
    channels_obj_one.getImageFileNames();
    channels_obj_one.sanityCheck();
    channels_obj_all = [channels_obj_all channels_obj_one];
    
    channels_obj_one = Channel([Analysis_Dir,filesep,all_rows(iM),'_Row_Analysis',filesep,'DAPI_Images']);
    channels_obj_one.getImageFileNames();
    channels_obj_one.sanityCheck();
    channels_obj_all = [channels_obj_all channels_obj_one];
    
    %% build MD
    
    % where to save it
    outputDirectory = [Analysis_Dir,filesep,all_rows(iM),'_Row_Analysis'];
    
    MD = MovieData(channels_obj_all,outputDirectory);
    MD.setPath(outputDirectory);
    MD.setFilename('movieData.mat')
    MD.sanityCheck();
    
    %% save movieData to movieList
    
    movie_cell{1,iMovie}= MD;
end

%% save movieList for the whole plate

mkdir([Analysis_Dir,filesep,'plate_movieList']);

ML = MovieList(movie_cell,[Analysis_Dir,filesep,'plate_movieList']);
ML.setPath([Analysis_Dir,filesep,'plate_movieList']);
ML.setFilename('movieList.mat')

ML.sanityCheck();

