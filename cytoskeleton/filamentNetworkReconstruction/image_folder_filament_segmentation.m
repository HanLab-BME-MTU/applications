function image_folder_filament_segmentation(image_folder, Parameter_MD, varargin)
% Function of single image filament segmentation with input MD from other
%               successfully segmented movie for the parameters in MD

% Input:      image_folder:    the folder of images to be segmented(one
%                              channel by one channel
%             Parameter_MD:    a loaded MD with good segmentation parameters
%                              if none, put [], so a default setting will be used
%                              with (1) Otsu with smoothing 1 (2) mask
%                              refine with 1 object, (3) image flatten with
%                              square, (4) steerable filter with [1 2], (5)
%                              segmentation with geo based alpha=2
%             pick_channel:    optional input, which channel to use in case there are more
%                               than one channel in the MD, default 1
%             keep_steps:      optional input, only ofr debugging, if 1,
%                               the steps are kept, default 0
%             output_dir:      optional input, if given, the results will
%                                be saved in this dir; if not given, will be in the image
%                                folder
%             whole_movie_filename: optional input, if given, the whole
%                                movie statistics will be loaded from this filename, which
%                                will be used in the filament segmentsation; by default,
%                                nothing is input, so whole movie stat will be calculated,
%                                well, in this case, without any meaningful effect

% Output:     this_MD: the movieData object, everything is saved as
%                           in GUI run

image_folder = GetFullPath(image_folder);
if(image_folder(end)~=filesep)
    image_folder = [image_folder,filesep];
end

tif_list = dir([image_folder, '*.tif']);
tiff_list = dir([image_folder, '*.tiff']);
jpg_list = dir([image_folder, '*.jpg']);

img_list = cell(1,1);
img_count = 0;
if numel(tif_list)>0
    for iT = 1 : numel(tif_list)
        img_count = img_count+1;
        img_list{img_count} = tif_list(iT).name;
    end
end

if numel(tiff_list)>0
    for iT = 1 : numel(tiff_list)
        img_count = img_count+1;
        img_list{img_count} = tiff_list(iT).name;
    end
end

if numel(jpg_list)>0
    for iT = 1 : numel(jpg_list)
        img_count = img_count+1;
        img_list{img_count} = jpg_list(iT).name;
    end
end

for iImage = 1 : img_count
    
    single_image_filament_segmentation(GetFullPath([image_folder,img_list{iImage}]), Parameter_MD, varargin{:});
    close all;
end
