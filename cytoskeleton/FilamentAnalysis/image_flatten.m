function movieData = image_flatten(movieData, paramsIn, varargin)

% for most cases, don't do background removal here.
background_removal_flag =0;

% Find the package of Filament Analysis
nPackage = length(movieData.packages_);

indexFilamentPackage = 0;
for i = 1 : nPackage
    if(strcmp(movieData.packages_{i}.getName,'FilamentAnalysis')==1)
        indexFilamentPackage = i;
        break;
    end
end

if(indexFilamentPackage==0)
    msgbox('Need to be in Filament Package for now.')
    return;
end


% Find the process of segmentation mask refinement.
nProcesses = length(movieData.processes_);

indexFlattenProcess = 0;
for i = 1 : nProcesses
    if(strcmp(movieData.processes_{i}.getName,'Image Flatten')==1)
        indexFlattenProcess = i;
        break;
    end
end

if indexFlattenProcess==0
    msgbox('Please set parameters for Image Flatten.')
    return;
end

% with no input funparam, use the one the process has on its own
if nargin < 2
    paramsIn = [];
    funParams = movieData.processes_{indexFlattenProcess}.funParams_;
else
    funParams = paramsIn;    
end

selected_channels = funParams.ChannelIndex;
flatten_method_ind = funParams.method_ind;
Gaussian_sigma = funParams.GaussFilterSigma;

TimeFilterSigma = funParams.TimeFilterSigma;
Sub_Sample_Num  = funParams.Sub_Sample_Num;

nFrame = movieData.nFrames_;

% default image flateen output dir
ImageFlattenProcessOutputDir = [movieData.outputDirectory_, filesep 'ImageFlatten'];

% if there is filamentanalysispackage
if (indexFilamentPackage>0)
    % and a directory is defined for this package
    if (~isempty(movieData.packages_{indexFilamentPackage}.outputDirectory_))
        % and this directory exists
        if (~exist(movieData.packages_{indexFilamentPackage}.outputDirectory_,'dir'))
            mkdir(movieData.packages_{indexFilamentPackage}.outputDirectory_);
        end
        ImageFlattenProcessOutputDir  = [movieData.packages_{indexFilamentPackage}.outputDirectory_, filesep 'ImageFlatten'];        
    end
end

if (~exist(ImageFlattenProcessOutputDir,'dir'))
    mkdir(ImageFlattenProcessOutputDir);
end

for iChannel = selected_channels
    ImageFlattenChannelOutputDir = [ImageFlattenProcessOutputDir,'/Channel',num2str(iChannel)];
    if (~exist(ImageFlattenChannelOutputDir,'dir'))
        mkdir(ImageFlattenChannelOutputDir);
    end
    
    output_dir_content = dir(fullfile([ImageFlattenChannelOutputDir,filesep,'*.*']));
    
    %if there are files in this dir, clear them
    if(length(output_dir_content)>2)
         delete([ImageFlattenChannelOutputDir,filesep,'*.*']);
    end
    
    movieData.processes_{indexFlattenProcess}.setOutImagePath(iChannel,ImageFlattenChannelOutputDir);
end

for iChannel = selected_channels
    for_temporal_filtering = cell(1);
    
    ImageFlattenProcessOutputDir = movieData.processes_{indexFlattenProcess}.outFilePaths_{iChannel};
    
    % Get frame number from the title of the image, this not neccesarily
    % the same as iFrame due to some shorting problem of the channel
    Channel_FilesNames = movieData.channels_(iChannel).getImageFileNames(1:movieData.nFrames_);
    
    filename_short_strs = uncommon_str_takeout(Channel_FilesNames);
    
    Frames_to_Seg = 1:Sub_Sample_Num:nFrame;
    Frames_results_correspondence = im2col(repmat(Frames_to_Seg, [Sub_Sample_Num,1]),[1 1]);
    Frames_results_correspondence = Frames_results_correspondence(1:nFrame);
    
    img_pixel_pool = [];
    for iFrame_subsample = 1 : length(Frames_to_Seg)
        iFrame = Frames_to_Seg(iFrame_subsample);
        currentImg = movieData.channels_(iChannel).loadImage(iFrame);
        
        if(size(currentImg,1)>100 || length(Frames_to_Seg)>20)
            smaller_currentImg = imresize(currentImg,[100 NaN]);
            if length(Frames_to_Seg)>50
                smaller_currentImg = imresize(currentImg,[50 NaN]);
            end
        else
            smaller_currentImg = currentImg;
        end
        
        img_pixel_pool = [img_pixel_pool smaller_currentImg(:)];
    end
    [hist_all_frame, hist_bin] = hist(double(img_pixel_pool),100);
    
    img_pixel_pool = double(img_pixel_pool(:));
    nonzero_img_pixel_pool= img_pixel_pool(img_pixel_pool>0);
    
    % if not found the loop use 3 times min
    low_005_percentile =  3*min(img_pixel_pool)+3*(max(img_pixel_pool)-min(img_pixel_pool))/100;
    for intensity_i = min(img_pixel_pool) : (max(img_pixel_pool)-min(img_pixel_pool))/100 : 3*min(img_pixel_pool)+3*(max(img_pixel_pool)-min(img_pixel_pool))/100
        if length(find(img_pixel_pool<=intensity_i))/length(img_pixel_pool)>0.005
            low_005_percentile = intensity_i;
            break;
        end
    end
    
    if(low_005_percentile==0 && flatten_method_ind==1)
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
    high_995_percentile = max(img_pixel_pool)/2;
    for intensity_i = max(img_pixel_pool) : -(max(img_pixel_pool)-min(img_pixel_pool))/100 : max(img_pixel_pool)/2
        if length(find(img_pixel_pool<intensity_i))/length(img_pixel_pool)<0.9995
            high_995_percentile = intensity_i;
            find_flag = 1;
            break;
        end
    end
    
    if exist('img_pixel_pool','var')
        delete 'img_pixel_pool';
    end
    
    %    low_005_percentile=0;
    %    high_995_percentile= 2^8-1;
    
    img_min=low_005_percentile;
    img_max=high_995_percentile;
    
       
    currentImg_cell = cell(1,1);
    
    for iFrame_subsample = 1 : length(Frames_to_Seg)
        hist_this_frame = hist_all_frame(:,iFrame_subsample);
        ind = find(hist_this_frame==max(hist_this_frame));
        center_value(iFrame_subsample) = hist_bin(ind(1));
        if(ind(1)>1)
            center_value_m1(iFrame_subsample) = hist_bin(ind(1)-1);
        else
            center_value_m1(iFrame_subsample)= center_value(iFrame_subsample);
        end
        %       center_value(iFrame_subsample) = 1;
    end
    center_value_int = mean((center_value_m1+center_value)/2);
    
    % record the stat numbers
    funParams.stat.low_005_percentile = low_005_percentile;
    funParams.stat.high_995_percentile = high_995_percentile;
    funParams.stat.center_value_int = center_value_int;
    
    
    center_value = center_value/max(center_value);
    center_value = sqrt(center_value);
    center_value = imfilter(center_value,[1 2 3 9 3 2 1]/21,'replicate','same');
    center_value = center_value/max(center_value);
    %     center_value(:)=1;
    img_min=0;
    img_max=high_995_percentile- center_value_int;
    
    % Make output directory for the flattened images
    ImageFlattenChannelOutputDir = movieData.processes_{indexFlattenProcess}.outFilePaths_{iChannel};
    if (~exist(ImageFlattenChannelOutputDir,'dir'))
        mkdir(ImageFlattenChannelOutputDir);
    end
    
    display('======================================');
    display(['Current movie: as in ',movieData.outputDirectory_]);
    display(['Start image flattening in Channel ',num2str(iChannel)]);
    
    for iFrame_subsample = 1 : length(Frames_to_Seg)
        iFrame = Frames_to_Seg(iFrame_subsample);
        disp(['Frame: ',num2str(iFrame)]);
        
        % Read in the intensity image.
        currentImg = movieData.channels_(iChannel).loadImage(iFrame);
        currentImg = double(currentImg);
        
        currentImg = currentImg - center_value_int;
        
        % Get rid of extreme noises
        currentImg(find(currentImg>high_995_percentile- center_value_int))=high_995_percentile- center_value_int;
        currentImg(find(currentImg<=0.00000001))=0.00000001;
        
        % based on the given method index, do log or sqrt to flatten the image
        if flatten_method_ind == 1
            
            currentImg = log(currentImg);
            currentImg = currentImg - log(img_min);
            currentImg = currentImg/...
                (log((center_value(iFrame_subsample))*img_max) ...
                -log((center_value(iFrame_subsample))*img_min));
            
        else
            
            currentImg = currentImg - img_min;
            currentImg = currentImg/(center_value(iFrame_subsample))/(img_max- img_min);
            if flatten_method_ind == 2
                currentImg = (currentImg).^(1/2);
            end
            
            if flatten_method_ind == 3
                currentImg = (currentImg).^(2/3);
            end
        end
        
        
        % Smooth the image in requested
        if Gaussian_sigma > 0
            currentImg = imfilter(currentImg, fspecial('gaussian',round(5*Gaussian_sigma), Gaussian_sigma),'replicate','same');
        end
       
        currentImg_cell{iFrame}=currentImg;
        
         %% %tif stack cost too much memory, geta cell instead
         if(TimeFilterSigma > 0)
             for_temporal_filtering{iFrame_subsample} = currentImg;
         else
             % if no need for time filtering, save to disk
             for sub_i = 1 : Sub_Sample_Num
                 if iFrame + sub_i-1 <= nFrame
                     
                     currentImg(currentImg<0)=0;
                     currentImg(currentImg>1)=1;                     
                     currentImg = currentImg*(2^8-1);
                     currentImg = uint8(currentImg);
                     
                     imwrite(currentImg,[ImageFlattenChannelOutputDir,'/flatten_', ...
                         filename_short_strs{iFrame + sub_i-1},'.tif']);
                 end
             end
         end
    end
    
    %% if temporal filtering is needed
    
    if(TimeFilterSigma > 0)
        
        % Initialize the after filtering cell
        after_temporal_filtering = cell(1);
        
        % prepare the temporal filter
        FilterHalfLength = 2*ceil(TimeFilterSigma);
        temperal_filter = zeros(1,1,2*FilterHalfLength+1);
        H = fspecial('gaussian',2*FilterHalfLength+1, TimeFilterSigma);
        H_1D = H(FilterHalfLength+1,:);
        H_1D = H_1D/(sum(H_1D));
        temperal_filter(1,1,:) = H_1D(:);
        
        % cut the sequence into trunks(length pace), to reduce the need of stacking a
        % long long sequence into one 3D matrix which kills the memory
        
        pace = min(length(Frames_to_Seg),10);
        
        for iTfilter = 1 : pace : length(Frames_to_Seg)
            
            % these are the index of starting and ending of image to be
            % filtered, so could include the neighboring image on both ends
            start_iF_sub = max(1, iTfilter-ceil(TimeFilterSigma));
            end_iF_sub = min(length(Frames_to_Seg), iTfilter+pace-1+ceil(TimeFilterSigma));
            
            % this is the results, without the padding
            start_filtered = iTfilter;
            end_filtered = min(length(Frames_to_Seg), iTfilter+pace-1);
            
            % initialize the 3D tensor
            image_tensor_10 = zeros(size(currentImg,1),size(currentImg,2),end_iF_sub-start_iF_sub+1);
            
            for iF_sub = start_iF_sub : end_iF_sub
                image_tensor_10(:,:,  iF_sub - start_iF_sub +1) ...
                    = for_temporal_filtering{iF_sub};
            end
            
            % filter the 3D tensor
            time_filtered_10 = imfilter(image_tensor_10,temperal_filter,'replicate','same');
            
            % read out into the results, without padding
            for iF_sub = start_filtered : end_filtered
                after_temporal_filtering{iF_sub}=...
                    squeeze(time_filtered_10(:,:,iF_sub-start_iF_sub+1));
            end
            
        end
        
        % save the filtered to hard disk
        for iFrame_subsample = 1 : length(Frames_to_Seg)
            iFrame = Frames_to_Seg(iFrame_subsample);
            disp(['Frame: ',num2str(iFrame)]);
            currentImg = after_temporal_filtering{iFrame_subsample};
            currentImg(currentImg<0)=0;
            currentImg(currentImg>1)=1;
            
            currentImg = currentImg*(2^8-1);           
            currentImg = uint8(currentImg);
            
            for sub_i = 1 : Sub_Sample_Num
                if iFrame + sub_i-1 <= nFrame
                    disp(['Frame: ',num2str(iFrame + sub_i-1)]);
                    
                    imwrite(currentImg,[ImageFlattenChannelOutputDir,'/flatten_', ...
                        filename_short_strs{iFrame + sub_i-1},'.tif']);
                end
            end
        end
    end
    
    %% if background removal is needed
    if background_removal_flag==1
        % Background substraction for uneven illumination
        for iFrame_subsample = 1 : length(Frames_to_Seg)
            iFrame = Frames_to_Seg(iFrame_subsample);
            disp(['Frame: ',num2str(iFrame)]);
            currentImg =  imread([ImageFlattenChannelOutputDir,'/flatten_', ...
                filename_short_strs{iFrame + sub_i-1},'.tif']);
            
            I = double(currentImg);
            % Get the x and y of the surface
            [XI, YI] = meshgrid(1:size(I,2), 1:size(I,1));
            % Fit a polynomial to the surface x-y both up to 2nd order
            fit_sur = fit([YI(:),XI(:)],I(:), 'poly22', 'Robust', 'on');
            % Reconstract the surface for uneven background---without the first part (the DC part)
            Z_fit = fit_sur.p01.*XI+ fit_sur.p10.*YI +fit_sur.p11.*XI.*YI+ ...
                fit_sur.p20.*YI.*YI+fit_sur.p02.*XI.*XI;
            % Substract the background
            currentImg = I-Z_fit;
            currentImg = currentImg/255;% back to 0~1 for saving image tif
            
            currentImg(currentImg<0)=0;
            currentImg(currentImg>1)=1;
            
            currentImg = currentImg*(2^8-1);          
            currentImg = uint8(currentImg);
            
            % Save to disk
            for sub_i = 1 : Sub_Sample_Num
                if iFrame + sub_i-1 <= nFrame
                    disp(['Frame: ',num2str(iFrame + sub_i-1)]);
                    
                    imwrite(currentImg,[ImageFlattenChannelOutputDir,'/flatten_', ...
                        filename_short_strs{iFrame + sub_i-1},'.tif']);
                end
            end
        end
    end
    
    %% 
    % this the end of "for" of each channel    
end
