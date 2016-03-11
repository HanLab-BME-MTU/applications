function [frames_points,frames_features] = get_MIP_points(path , file_name, MIP_params)

%% ================================================================================
% get_MIP_points 
% -------------------------------------------------
% Get the MIP feature points vector for this clip
%
% Inputs:
% ~~~~~~
% path - full path to the video file
% file_name - the name of the video file
% MIP_params - The parameters for MIP coding  (params.MIP_params)
%
% Outputs:
% ~~~~~~~~
% frames_points - matrix points [x,y,t] of sampled frames
% frames_features - matrix of histograms of sampled frames (histogram per patch)
%
%
% Copyright 2012, Orit Kliper-Gross, Yaron Gurovich, Tal Hassner, and Lior Wolf
%
%% ================================================================================

% Read movie:
try
   movObj = mmreader( fullfile(path,file_name) );
   nFrames = movObj.NumberOfFrames;
catch m
    disp(m);
    error(['Cannot open file ',fullfile(path,file_name)]);
end
    
% Init out:
frames_features = {};   % histograms per sampled frame
frames_points = {}; % points [x,y,t] per sampled frame
    
% define the frames id to process:
frames_to_process = 1+MIP_params.movment_int : MIP_params.frame_gap : nFrames-MIP_params.movment_int;

% Read first frame for initializing:
Current_F = read(movObj,1);
if(MIP_params.resize_image == 1) 
    % Set scalling:
    if(isfield(MIP_params,'resize_image_value') && (MIP_params.resize_image_value > 0))
        rescale = MIP_params.resize_image_value;
    else
        rescale = 100 / size(Current_F,1);
    end
    % gray only:
    if size(Current_F,3)>1
        Current_F = rgb2gray(Current_F);
    end
    Current_F = single(Current_F);
    if (rescale < 0.8) || (isfield(MIP_params,'resize_image_value') && (MIP_params.resize_image_value > 0))
        Current_F = imresize(Current_F, rescale);
    else
        rescale = 1;
    end
end
    
% Set val_map:(for coding bits per pixel)
val_map = zeros(size(Current_F,1),size(Current_F,2),8);
bit_val = [1 2 4 8 16 32 64 128]';
for i = 1: size(val_map,1)
    for j = 1: size(val_map,2)
        val_map(i,j,:) = bit_val;
    end
end

% Debug mode preperation
% ---------------------------------

% Debug (coding) mode preperation:
if (MIP_params.MIP_dbg_dumps)
    % debug into mov name directory:
    if(isempty(MIP_params.dbg_dir))
        MIP_params.dbg_dir = fullfile('..','..','DebugDir');            
    end
    MIP_params.dbg_dir = fullfile(MIP_params.dbg_dir,file_name(1:end-4));
    if(~exist(MIP_params.dbg_dir,'dir'))
        mkdir(MIP_params.dbg_dir);
    end          
    tmp_frame = fullfile(MIP_params.dbg_dir,['dbg_cur_frame','.jpg']);
    % write avi with Xvid compression
    o_name = fullfile(MIP_params.dbg_dir,['dbg_',file_name(1:end-4),'.avi']);
    if(~MIP_params.MIP_dbg_dumps_imonly)
        warning off;
        hvfw = avifile(o_name,'fps',movObj.FrameRate,'compression','Xvid');
    end
end

% Debug (align) mode preperation:
if(MIP_params.align.en && MIP_params.align.dbg_out)                
    if(isempty(MIP_params.align.dbg_dir)) % should not be empty (filled in save_MIP_features_directions_coding)
        MIP_params.align.dbg_dir = fullfile('..','..','Data','DebugDir','debug_align');            
    end
    if(~exist(MIP_params.align.dbg_dir,'dir'))
        mkdir(MIP_params.align.dbg_dir);
    end              
    o_name = fullfile(MIP_params.align.dbg_dir,['MIP_Align_',file_name(1:end-4),'.avi']);        
    warning off;
    hvfwS = avifile(o_name,'fps',movObj.FrameRate,'compression','Xvid');
end
       
% Precomputed wrapped frames support:
% ---------------------------------------------------

load_precomp_wrap_data = false;
if( MIP_params.align.en && MIP_params.align.precompute.en )
    if(isempty(MIP_params.align.precompute.dir_name)) % should not be empty (filled in save_MIP_features_directions_coding)
        MIP_params.align.precompute.dir_name = fullfile('..','..','Data','OutputDir','wrapped_frames');% If no location store in the default OutputDir
    end
    precomp_full_name = fullfile(MIP_params.align.precompute.dir_name ,[file_name(1:end-4),'_wrapped.mat']);
    % load existing file if no request to override
    if(exist(precomp_full_name,'file') && ~MIP_params.align.precompute.overide)
        precomp_wrap_data = load(precomp_full_name);
        precomp_wrap_data = precomp_wrap_data.precomp_wrap_data;
        load_precomp_wrap_data = true;
    end
end

% main loop (on relevant frames):
% -----------------------------------------
MIP_params.valid_mask_en = false;

for f_id = 1:length(frames_to_process)
    
    % read frames:
    f = frames_to_process(f_id);                
    Prev_F    = read(movObj,f - MIP_params.movment_int);% mov(f - MIP_params.movment_int).cdata;%
    Current_F = read(movObj,f);% mov(f).cdata;%
    Next_F    = read(movObj,f + MIP_params.movment_int);% mov(f + MIP_params.movment_int).cdata;%

    % high res images:
    Prev_F_org = Prev_F;
    Current_F_org = Current_F;
    Next_F_org = Next_F;

    if size(Current_F,3)>1 % if colored movie change to gray levels
        Prev_F = rgb2gray(Prev_F);
        Current_F = rgb2gray(Current_F);
        Next_F = rgb2gray(Next_F);
    end
    Prev_F = single(Prev_F);
    Current_F = single(Current_F);
    Next_F = single(Next_F);

    % resize image:
    if((MIP_params.resize_image == 1) && (rescale ~= 1))
        Prev_F = imresize(Prev_F, rescale);
        Current_F = imresize(Current_F, rescale);
        Next_F = imresize(Next_F, rescale);
    end

    % alignment:
    % --------------
    
    if(MIP_params.align.en)
        % If no stored wrapping, compute it:
        if(~MIP_params.align.precompute.en || ~load_precomp_wrap_data)
            [wrap_data] = frame_MIP_align(Prev_F, Current_F, Next_F,Prev_F_org, Current_F_org, Next_F_org,rescale,MIP_params);
        end
        % load/store warpping:
        if(MIP_params.align.precompute.en)
            if(load_precomp_wrap_data)%load
                wrap_data = precomp_wrap_data{f_id};
                if(isempty(wrap_data))
                    error('Not consintent data , cant find relevant frame');
                end
            else % store
                precomp_wrap_data{f_id} = wrap_data;
            end
        end
        Next_F = wrap_data.Next_F_w;
        Prev_F = wrap_data.Prev_F_w;
        MIP_params.valid_mask = wrap_data.valid_mask;
        MIP_params.valid_mask_en = true;    
        
        % Debug (alignemnt) out mov
        % -----------------------------------
        if(MIP_params.align.dbg_out)% Wrap Org mov and store to a video stream:
            Next_F_w = double(wrap_data.Next_F_w);
            if (rescale < 1)% get bigger frame for XVID compatibility
                Next_F_w = imresize(Next_F_w,1/rescale);
            end
            I = zeros([size(Next_F_w,1),size(Next_F_w,2),3],'double');
            for in=1:3
                I(:,:,in) = Next_F_w;
            end
            I = im2uint8(I/max(I(:)));
            if(min(size(I)) < 480)%for xvid 
                I = imresize(I,[480,640]);
            end
            hvfwS = addframe(hvfwS,I);
        end
    end % if (MIP_params.align.en)
    
    % Main features computation for the current frame:
    % ----------------------------------------------------------------
    [frames_points{f_id},frames_features{f_id},MIP_for_display1, MIP_for_display2] = ...
        frame_MIP_features(f,Prev_F, Current_F, Next_F, val_map, MIP_params); 

    % Debug (coding)
    % --------------------
    if(MIP_params.MIP_dbg_dumps)
        % Prepare coded frames to plot:
        Cur = (((Current_F/max(Current_F(:)))*255));
        disp_img = zeros( size(Current_F,1), size(Current_F,2), 3, 'single');
        disp_img(:,:,1) = floor(Cur + MIP_for_display1 * (255/8) - MIP_for_display2 * (255/8));
        disp_img(:,:,2) = floor(Cur - MIP_for_display1 * (255/8) - MIP_for_display2 * (255/8));
        disp_img(:,:,3) = floor(Cur - MIP_for_display1 * (255/8) + MIP_for_display2 * (255/8));
        h = figure('visible','off');
        disp_img = im2uint8(disp_img/max(disp_img(:)));
        image(disp_img);
        title(['MIP of frame #' num2str(f),' points # ',num2str(size(frames_points{f_id},2))]);
        axis off;
        saveas(h,[tmp_frame(1:end-4),'_',num2str(f),'.jpg'],'jpg');
        close(h);            
        if (rescale < 1)% get bigger frame for XVID compatibility
            disp_img = imresize(disp_img,1/rescale);
        end
        if(~MIP_params.MIP_dbg_dumps_imonly)
            if(min(size(disp_img)) < 480)%for xvid 
                disp_img = imresize(disp_img,[480,640]);
            end

            hvfw = addframe(hvfw,disp_img);
        end
    end % if(MIP_params.MIP_dbg_dumps) ...

end % for f_id...

% Close debug mov:
if (MIP_params.MIP_dbg_dumps && ~MIP_params.MIP_dbg_dumps_imonly)
    hvfw = close(hvfw);
end
if(MIP_params.align.en) && (MIP_params.align.dbg_out)
    hvfwS = close(hvfwS);        
end

% Store wrap frames:
if(MIP_params.align.en && MIP_params.align.precompute.en && ~load_precomp_wrap_data)
    if(~exist(MIP_params.align.precompute.dir_name ,'dir'))
        mkdir(MIP_params.align.precompute.dir_name);
    end
    save(precomp_full_name,'precomp_wrap_data','-v7.3');
end

% Convert features to a one matrix for all frames:
frames_points = (cat(2,frames_points{:}))';
frames_features = (cat(2,frames_features{:}))';

end % of function: get_MIP_points