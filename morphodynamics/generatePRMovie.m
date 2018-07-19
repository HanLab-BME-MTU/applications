function generatePRMovie(handles)
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

if ~isfield(handles, 'snapShot')
    snapShot = 0;
else
    snapShot = handles.snapShot;
end

if snapShot == 1
    if ~exist([output_dir filesep 'snapShot'], 'dir')
       mkdir(output_dir, 'snapShot');
    end
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

mask1 = imread([handles.directory_name filesep maskdir(1).name]);
[rowSize, colSize] = size(mask1);
MAP_SCALE = min(round(rowSize/2),100);
cmap_edge_evolution=jet(MAP_SCALE);

nColumn = handles.maskinfo.mask_num-1;

QTMovieName{1} = [output_dir filesep 'Activity_movie.mov'];
if CC > 1 && (ISCLOSE  == 1)
    QTMovieName{2} = [output_dir filesep 'Activity_movie_CC2_RES.mov'];
end





for q_idx = 1:length(QTMovieName)
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
            row_diff = centroid(time_idx,:) - centroid(1,:);
            center_diff = repmat(row_diff, [size(this_pixel_edge{time_idx},1) 1]);
            this_pixel_edge{time_idx} = this_pixel_edge{time_idx} + center_diff;
        end        
    end
    
    if (BESTSTART == 1) && (ISCLOSE == 1)
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

    MAX_COLOR_VAL_PROT = abs(robust_av_protrusion(1)+0.5*robust_av_protrusion(2)*prot_x_time_axis(end)) + 3*robust_sigma_protrusion;
    protrusion_normal(protrusion_normal > MAX_COLOR_VAL_PROT) = MAX_COLOR_VAL_PROT;
    protrusion_normal(protrusion_normal <-MAX_COLOR_VAL_PROT) = -MAX_COLOR_VAL_PROT;
    MAP_COLOR = round(((protrusion_normal + MAX_COLOR_VAL_PROT)/(2*MAX_COLOR_VAL_PROT)) * (MAP_SCALE-1))+ 1;

    scalebar = zeros(MAP_SCALE,10,3);
    for i=1:MAP_SCALE
        scalebar(i,1:10,1) = cmap_edge_evolution(MAP_SCALE-i+1,1);
        scalebar(i,1:10,2) = cmap_edge_evolution(MAP_SCALE-i+1,2);
        scalebar(i,1:10,3) = cmap_edge_evolution(MAP_SCALE-i+1,3);
    end

    if isfield(handles, 'imgDir')
        imgDir = dir([handles.imgDir filesep '*.tif']);
    end

    MakeQTMovie('start', QTMovieName{q_idx});

    for time_idx = 1:nColumn
        pixel_list = this_pixel_edge{time_idx};
        round_pixel_list = round(pixel_list);
        pixelNum = size(round_pixel_list,1);

        nfold = SEG_NR;
        binsize = floor(pixelNum/nfold);
        start_idx = 0;
        idx{1} = mod(start_idx:(start_idx+binsize-1), pixelNum)+1;
        start_idx = mod(start_idx+binsize, pixelNum)+1;
        for w=2:(nfold-mod(pixelNum, nfold))
            idx{w} = start_idx:(start_idx+binsize-1);
            start_idx  = start_idx+binsize;
        end
        for w=(nfold+1-mod(pixelNum, nfold)):nfold
            idx{w} = start_idx:(start_idx+binsize);
            start_idx  = start_idx+binsize+1;
        end
        sector_idx = [];
        for j=1:nfold
            sector_idx = [sector_idx; j*ones(length(idx{j}),1)];
        end

        if isfield(handles, 'imgDir')
            img_edge_rgb = imread([handles.imgDir filesep imgDir(time_idx).name]);
            img_edge_rgb = repmat(mat2gray(img_edge_rgb), [1,1,3]);
        else
            img_edge_rgb = zeros(rowSize,colSize,3);
        end
        for j=1:pixelNum
            if round_pixel_list(j,1) > 0 && round_pixel_list(j,1) <= colSize && round_pixel_list(j,2) > 0 && round_pixel_list(j,2) <= rowSize
                img_edge_rgb(round_pixel_list(j,2),round_pixel_list(j,1),1)=cmap_edge_evolution(MAP_COLOR(sector_idx(j),time_idx),1);
                img_edge_rgb(round_pixel_list(j,2),round_pixel_list(j,1),2)=cmap_edge_evolution(MAP_COLOR(sector_idx(j),time_idx),2);
                img_edge_rgb(round_pixel_list(j,2),round_pixel_list(j,1),3)=cmap_edge_evolution(MAP_COLOR(sector_idx(j),time_idx),3);

                img_edge_rgb(1+25:MAP_SCALE+25,[1:10]+10,:) = scalebar;
            end
        end

        h_QT = figure(101);
        imshow(img_edge_rgb, []);
        title('Protrusion/Retraction in Sectors');
        text(5,10, ['#' num2str(time_idx)], 'Color', 'white', 'FontSize', 12);
        text(22,25, num2str(round(MAX_COLOR_VAL_PROT*10)/10), 'Color', 'yellow', 'FontSize', 9);
        text(22,25 + MAP_SCALE, num2str(-round(MAX_COLOR_VAL_PROT*10)/10), 'Color', 'yellow', 'FontSize', 9);
        hold on;
        plot(round_pixel_list(1,1), round_pixel_list(1,2),'ro','MarkerFaceColor','r','MarkerSize',5);
        hold off;

        MakeQTMovie addfigure
        if snapShot
            print(h_QT, [output_dir filesep 'snapShot' filesep 'snap' num2str(time_idx) '.eps'],'-depsc2');
        end
    end
    close(h_QT);

    MakeQTMovie('finish');
    fprintf(1,'%s has been generated.\n',QTMovieName{q_idx});
end
