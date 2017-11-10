function image_flatten_folder(input_folder, flatten_method,varargin)
% Function for flatten images ( in the folder setting, not the movieData setting)
% Input: input_folder: the folder consisting the images to be flattened.
%        flatten_method: a string for the method of flatten: 
%                       'log': taking the log
%                       'sqrt': taking the sqrt( if anything other than 'log' is inputed, 'sqrt' will be used.)
%                       note all output images are normalized to 0~1, with extreme intensity pixel correction
%        Optional input: output_folder, if not given, will be named flatten_image side by side with the input folder
% No output, images will be saved to the disk.

% Liya Ding 2013.

% parse input
ip = inputParser;
ip.addRequired('input_folder', @ischar); 
ip.addRequired('flatten_method', @ischar);
ip.addOptional('output_folder', 'flatten_image', @ischar);
ip.parse(input_folder, flatten_method, varargin{:});
output_folder = ip.Results.output_folder;

%set location and make output dir
outputLocation = fileparts(input_folder);

if(isempty(outputLocation))
    outputLocation='.';
end

outputDir = ([outputLocation filesep output_folder]);

if(~exist(outputDir,'dir'))
    mkdir(outputDir);
end

%read in files and show in window
[input_files,~,input_index] = imDir(input_folder);
disp('input_files = ' );
disp(input_files);

nFrame = length(input_index);

img_pixel_pool = [];
for iFrame = 1 : nFrame
    currentImg = double(imread([input_folder filesep input_files(iFrame).name]));
    img_pixel_pool = [img_pixel_pool currentImg(:)];
end

img_pixel_pool = double(img_pixel_pool(:));

nonzero_img_pixel_pool = img_pixel_pool(img_pixel_pool~=0);

% if not found the loop use 3 times min
low_005_percentile =  3*min(img_pixel_pool)+3*(max(img_pixel_pool)-min(img_pixel_pool))/100;
for intensity_i = min(img_pixel_pool) : (max(img_pixel_pool)-min(img_pixel_pool))/100 : 3*min(img_pixel_pool)+3*(max(img_pixel_pool)-min(img_pixel_pool))/100
    if length(find(img_pixel_pool<=intensity_i))/length(img_pixel_pool)>0.005
        low_005_percentile = intensity_i;
        break;
    end
end

if(low_005_percentile==0 && strcmp(flatten_method,'log'))
    for intensity_i = min(img_pixel_pool) : (max(img_pixel_pool)-min(img_pixel_pool))/100 : 3*min(img_pixel_pool)+3*(max(img_pixel_pool)-min(img_pixel_pool))/100;
        if length(find(img_pixel_pool<=intensity_i))/length(img_pixel_pool)>0.01
            low_005_percentile = intensity_i;
            break;
        end
    end
    if(low_005_percentile==0)
        low_005_percentile = 3*min(nonzero_img_pixel_pool);
    end
end

% if not found the loop use half max
high_9995_percentile = max(img_pixel_pool)/2;
for intensity_i = max(img_pixel_pool) : -(max(img_pixel_pool)-min(img_pixel_pool))/100 : max(img_pixel_pool)/2
    if length(find(img_pixel_pool<intensity_i))/length(img_pixel_pool)<0.9995
        high_9995_percentile = intensity_i;
        find_flag = 1;
        break;
    end
end

img_min=low_005_percentile;
img_max=high_9995_percentile;

[hist_all_frame, hist_bin] = hist(img_pixel_pool,5);
    
    display('Start image flattening in the folder');
    
    for iFrame = 1 : nFrame
        disp(['Frame: ',num2str(iFrame)]);
        
        % Read in the intensity image.
        currentImg = double(imread([input_folder filesep input_files(iFrame).name]));
        % Get rid of extreme noises
        currentImg(find(currentImg<low_005_percentile))=low_005_percentile;
        currentImg(find(currentImg>high_9995_percentile))=high_9995_percentile;
        
        % based on the given method index, do log or sqrt to flatten the image
        if(strcmp(flatten_method,'log'))
            % taking the log on 
            currentImg = log(currentImg);
            currentImg = currentImg - log(img_min);
            currentImg = currentImg/(log(img_max)-log(img_min));
        else
            currentImg = currentImg - img_min;
            currentImg = currentImg/(img_max- img_min);
            currentImg = (currentImg).^(1/2);
        end
        
        imwrite(currentImg,[outputDir filesep input_files(iFrame).name '_flatten.tif']);
        
    end
end
