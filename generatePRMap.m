function generatePRMap(handles)
% generatePRMap generates the protrusion Activity map
%
% Input:    handles.batch_processing        flag for batch processing
%           handles.directory_name          mask directory
%           handles.FileType                mask format
%           handles.result_directory_name   result directory
%           handles.timevalue               time resolution
%           handles.resolutionvalue         spatial resolution
%           handles.segvalue                number of segment for activity map
%           handles.dl_rate                 down sample to # of pixels if length for fragment is greater than #
%           handles.maskinfo                all the information about mask
%
% See also: samPrPanel, protrusionAnalysis
% Last updated: August 26, 2008 by Shann-Ching Chen, LCCB

if  ~isfield(handles,'CC')
    % centroid correction
    CC = 0;
else
    CC = handles.CC;
end

if ~isfield(handles,'BESTSTART')
    BESTSTART = 0;
else
    BESTSTART = handles.BESTSTART;
end

if CC > 1
    handles.outdir = [handles.result_directory_name  filesep 'analysis_CC2_dl' num2str(handles.dl_rate)];
else
    handles.outdir = [handles.result_directory_name  filesep 'analysis_dl' num2str(handles.dl_rate)];
end

output_dir = handles.outdir;

if ~isfield(handles, 'FileType')
    fullmaskpath = [handles.directory_name filesep '*.tif'];
else
    fullmaskpath = [handles.directory_name filesep handles.FileType];
end

if ~isfield(handles, 'ISCLOSE')
    ISCLOSE = 0;
else
    ISCLOSE = handles.ISCLOSE;
end

% file name of the the pixel edge
file_pixel_edge = [output_dir filesep 'pixel_edge.mat'];

% file name of the normal vectors
file_normal=[output_dir filesep 'normal_matrix.mat'];

% get the filelist filenames containing the protrusion vectors
file_protrusion=[output_dir filesep 'protrusion.mat'];

% get the filelist filenames containing the centroid coordinates
file_centroid=[output_dir filesep 'centroid.mat'];

%%%%%%%%%% load data  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load(file_pixel_edge);
load(file_normal);
load(file_protrusion);
load(file_centroid);

if BESTSTART == 1
    output_dir = [handles.outdir 'B1'];
elseif BESTSTART == -1
    output_dir = [handles.outdir 'B-1'];
end

if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end

% folder for the result figures
if ~exist([output_dir filesep 'figures'], 'dir')
    mkdir(output_dir, 'figures');
end

% folder for the result figures
if ~exist([output_dir filesep 'result_mat'], 'dir')
    % folder for the result .mat
    mkdir(output_dir, 'result_mat');
end

% number of segment
SEG_NR        = handles.segvalue;
PIXEL         = handles.resolutionvalue;
TIME_INTERVAL = handles.timevalue;
UNITS = 0;
if UNITS == 0
    xLabelText = 'frame #';
    yLabelText = 'Velocity (pixel/frame)';
elseif UNITS == 1
    xLabelText = 'Time (s)';
    yLabelText = 'Velocity (nm/s)';
else
    xLabelText = 'Time (min)';
    yLabelText = 'Velocity (um/min)';
end


maskdir = dir(fullmaskpath);
mask_num = length(maskdir);
maskinfo.fullmaskpath = fullmaskpath;
maskinfo.maskdir = maskdir;
maskinfo.mask_num = mask_num;
handles.maskinfo = maskinfo;






nColumn = handles.maskinfo.mask_num-1;

img_prot_filename{1} = 'img_prot.mat';
prot_activity_tif_filename{1} = 'prot_activity_img.tif';
prot_activity_eps_filename{1} = 'prot_activity_img.eps';
if CC > 1 & (ISCLOSE  == 1)
    img_prot_filename{2} = 'img_protCC2.mat';
    prot_activity_tif_filename{2} = 'prot_activity_imgCC2.tif';
    prot_activity_eps_filename{2} = 'prot_activity_imgCC2.eps';    
end

for q_idx = 1:length(img_prot_filename)
    protrusion_normal = zeros(SEG_NR, nColumn);
    if q_idx == 1
        this_protrusion = protrusion;        
        this_normal_matrix = normal_matrix;
        this_pixel_edge = pixel_edge;        
    else
        this_protrusion = protrusion_CC2;
        this_normal_matrix = normal_matrix;
        this_pixel_edge = pixel_edge;        
        for time_idx = 2:length(this_protrusion)
            row_diff = [centroid(time_idx,:) - centroid(1,:)];
            center_diff = repmat(row_diff, [size(this_pixel_edge{time_idx},1) 1]);
            this_pixel_edge{time_idx} = this_pixel_edge{time_idx} + center_diff;
        end   
    end
    
    if (BESTSTART == 1) & (ISCLOSE == 1)
        for pixel_idx = 1:size(this_protrusion{1},1)
            startPixel_InFrame = [];
            for time_idx = 1:mask_num-1
                if pixel_idx == 1
                    % good for starting position
                    pT = this_protrusion{time_idx}(pixel_idx,3:4);
                    nT = this_normal_matrix{time_idx}(pixel_idx,3:4);
                    tmpI = 1;
                else
                    % for other positions, find correspondence
                    if time_idx == 1
                        startposTM1 = this_protrusion{1}(pixel_idx,1:2);
                        startposT   = startposTM1 + this_protrusion{1}(pixel_idx,3:4);
                        pT = this_protrusion{time_idx}(pixel_idx,3:4);
                        nT = this_normal_matrix{time_idx}(pixel_idx,3:4);
                        %                     VM1 = startposTM1 - centroid(1,:);
                        %                     V = startposT - centroid(2,:);
                        tmpI = pixel_idx;
                    else
                        % In protrusion{time_idx}, what is the closest node to startposT?
                        [tmp, tmpI] = min(sum((repmat(startposT, [size(this_protrusion{time_idx},1) 1]) - this_protrusion{time_idx}(:,1:2)).^2,2));
                        pT = this_protrusion{time_idx}(tmpI,3:4);
                        nT = this_normal_matrix{time_idx}(tmpI,3:4);
                        startposTM1 = this_protrusion{time_idx}(tmpI,1:2);
                        startposT = this_protrusion{time_idx}(tmpI,1:2) + this_protrusion{time_idx}(tmpI,3:4);
                    end
                end
                startPixel_InFrame = [startPixel_InFrame tmpI];

                nT = nT./norm(nT);
                scal            = pT(1) * nT(1) + pT(2) * nT(2);
                x_av_prot_proj  = scal .* nT(1);
                y_av_prot_proj  = scal .* nT(2);
                prot_proj(time_idx) = sign(scal).*sqrt(x_av_prot_proj.^2 + y_av_prot_proj.^2);
            end
            startPixelArray{pixel_idx} = startPixel_InFrame;
            protArray{pixel_idx} = prot_proj;
            meanProtArray(pixel_idx) = mean(prot_proj);
        end
        [minY, minI] = min(meanProtArray);
        for time_idx = 1:size(startPixelArray{pixel_idx},2)
            this_protrusion{time_idx} = [this_protrusion{time_idx}(startPixelArray{minI}(time_idx):end,:);...
                this_protrusion{time_idx}(1:startPixelArray{minI}(time_idx)-1,:)];

            this_normal_matrix{time_idx} = [this_normal_matrix{time_idx}(startPixelArray{minI}(time_idx):end,:);...
                this_normal_matrix{time_idx}(1:startPixelArray{minI}(time_idx)-1,:)];

            this_pixel_edge{time_idx} = [this_pixel_edge{time_idx}(startPixelArray{minI}(time_idx):end,:);...
                this_pixel_edge{time_idx}(1:startPixelArray{minI}(time_idx)-1,:)];
        end
    end

    for time_idx = 1:nColumn
        % average normal vectors over number of segment
        [x_normal_out_av, y_normal_out_av, x_av_pos_normal, y_av_pos_normal]=...
            prGetAvEdgeSam(this_normal_matrix{time_idx}, 'nr_sect', SEG_NR);

        % re-normalize the averaged normals
        l_n_av = sqrt(x_normal_out_av.^2 + y_normal_out_av.^2);
        x_av_normal = x_normal_out_av ./ l_n_av;
        y_av_normal = y_normal_out_av ./ l_n_av;

        % average protrusion vectors over number of segment
        [x_av_prot, y_av_prot, x_av_pos_prot, y_av_pos_prot]=...
            prGetAvEdgeSam(this_protrusion{time_idx}, 'nr_sect', SEG_NR);

        av_prot         = sqrt(x_av_prot.^2 + y_av_prot.^2);
        scal            = x_av_normal .* x_av_prot + y_av_normal .* y_av_prot;
        x_av_prot_proj  = scal .* x_av_normal;
        y_av_prot_proj  = scal .* y_av_normal;
        av_prot_proj    = sign(scal).*sqrt(x_av_prot_proj.^2 + y_av_prot_proj.^2);
        protrusion_normal(:,time_idx) = av_prot_proj' * PIXEL / TIME_INTERVAL;
    end

    protTimePoints = [1:nColumn];
    prot_x_time_axis = protTimePoints.* TIME_INTERVAL;
    av_protrusion_n_time = sum(protrusion_normal, 1) / SEG_NR;
    [robust_av_protrusion, stats] = robustfit(prot_x_time_axis, av_protrusion_n_time);
    robust_sigma_protrusion = stats.robust_s;
%%%%%%%%%%%%%%%%  From here, it's different from eneratePRMovie.m
    IMAGE_STRETCH   = 15;
    IMG_GAUSS_SIG   = 6;
    project_img_protrusion = protrusion_normal;
    img_protrusion  = imresize(project_img_protrusion,IMAGE_STRETCH, 'nearest');
    [img_protrusion_f, dum] = Gauss2DBorder(img_protrusion, IMG_GAUSS_SIG);

    h_protrusion_img = figure('Visible', 'off');
    imagesc(img_protrusion_f);
    V = axis;
    No_Ticks = 10;
    if No_Ticks > nColumn
        No_Ticks = nColumn;
    end

    XTickV_STEP = (V(1,2)-V(1,1))/No_Ticks;
    XTickV = V(1,1):XTickV_STEP:V(1,2);
    YTickV_STEP = (V(1,4)-V(1,3))/No_Ticks;
    YTickV = V(1,3):YTickV_STEP:V(1,4);
    XTickL = {};
    YTickL = {};
    for jj = 1:length(XTickV)
        XTickL = [XTickL num2str(floor(XTickV(jj)/IMAGE_STRETCH)*TIME_INTERVAL)];
    end
    for jj = 1:length(YTickV)
        YTickL = [YTickL num2str(floor(YTickV(jj)/IMAGE_STRETCH))];
    end
    set(gca,'XTick', XTickV);
    set(gca,'XTickLabel',XTickL);
    set(gca,'YTick', YTickV);
    set(gca,'YTickLabel',YTickL);
    xlabel(xLabelText);
    ylabel('Sector #');
    title(['Protrusion ' yLabelText]);
    colormap(jet);
    MAX_COLOR_VAL_PROT = abs(robust_av_protrusion(1)+0.5*robust_av_protrusion(2)*prot_x_time_axis(end)) + 3*robust_sigma_protrusion;
    caxis([-MAX_COLOR_VAL_PROT MAX_COLOR_VAL_PROT]);
    colorbar;
    
    save([output_dir filesep 'result_mat' filesep img_prot_filename{q_idx}], 'img_protrusion', 'img_protrusion_f', 'project_img_protrusion','MAX_COLOR_VAL_PROT');
    mapPath = [output_dir filesep 'figures' filesep prot_activity_tif_filename{q_idx}];
    epsPath = [output_dir filesep 'figures' filesep prot_activity_eps_filename{q_idx}];
    print(h_protrusion_img, mapPath,'-dtiff');
    print(h_protrusion_img, epsPath,'-depsc2');

    if handles.batch_processing
        % added for batch processing, Sam, 08/22/2008
        fprintf(1,'%s\n',['Activity Map is located at ' mapPath]);
    else
        set(handles.statusbar, 'String', ['Activity Map is located at ' mapPath '.']);
    end
end