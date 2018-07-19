function [OK, handles] = protrusionAnalysis(handles)
% protusionAnalysis measures the protrusion rate of cell leading edge based
% on constraint optimization.
%
% Input:    handles.batch_processing        flag for batch processing
%           handles.directory_name          mask directory
%           handles.FileType                mask format
%           handles.result_directory_name   result directory
%           handles.timevalue               time resolution, seconds per frame
%           handles.resolutionvalue         spatial resolution, nms per pixel
%           handles.segvalue                number of segments for activity map
%           handles.dl_rate                 down sample to # of pixels if length for fragment is greater than #
%
% Output:   handles.maskinfo                all the information about mask
%           handles                         other fields have the same structure as input fields
%
% Last updated: Octobor 08, 2008 by Shann-Ching Chen, LCCB, with CC and ratio added
% See also: pandaPanel, generatePRMap
% handles.ISCLOSE = 1;
% generatePRMovie(handles);
% keyboard;

OK = 1;
ISCLOSE = 1;

if handles.batch_processing > 0;   % added for batch processing, Sam, 08/22/2008
    % added for central correction, Sam, 09/02/2008
    % when batch_processing == 0, no batch
    % when batch_processing == 1, batch, no CC
    % when batch_processing == 2, batch, CC
    CC = handles.batch_processing;
    fullmaskpath = [handles.directory_name filesep handles.FileType];
    maskdir = dir(fullmaskpath);
    mask_num = length(maskdir);
    maskinfo.fullmaskpath = fullmaskpath;
    maskinfo.maskdir = maskdir;
    maskinfo.mask_num = mask_num;
    maskinfo.FileType = handles.FileType;
    handles.maskinfo = maskinfo;
    firstmask = imread([handles.directory_name filesep handles.maskinfo.maskdir(1).name]);
    handles.rowSize = size(firstmask,1);
    handles.colSize = size(firstmask,2);
    handles.CC = CC;
else
    if ~isfield(handles,'CC')
        CC = 1;
    else
        CC = handles.CC;
    end
end

if ~isfield(handles,'ORIENT_CELL')
    %smoothing parameter of the approximating spline
    ORIENT_CELL = 4;
else
    ORIENT_CELL = handles.ORIENT_CELL;
end

% in pandaParameter.m
% if get(handles.optimized,     'Value');
%     parameters.BESTSTART  = 1;
% elseif get(handles.follow,    'Value');
%     parameters.BESTSTART  = 0;
% elseif get(handles.fixed,     'Value');
%     parameters.BESTSTART  = -1;
% else
%     parameters.BESTSTART  = 0;
% end

if ~isfield(handles,'BESTSTART')
    % default: starting pixel follows protrusion vector
    handles.BESTSTART = 0;
end

%%%%%% Parameters for the spline function  %%%%%%
if ~isfield(handles,'TOLERANCE')
    %smoothing parameter of the approximating spline
    TOLERANCE=30;
else
    TOLERANCE = handles.TOLERANCE;
end

if CC > 1
    % added for central correction, Sam, 09/02/2008
    handles.outdir = [handles.result_directory_name  filesep 'analysis_CC2_dl' num2str(handles.dl_rate)];
else
    handles.outdir = [handles.result_directory_name  filesep 'analysis_dl' num2str(handles.dl_rate)];
end
output_dir = handles.outdir;

if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end

exist_flag = exist(output_dir,'dir');
if exist_flag == 0
    %it does not exist, so try to create it
    [s, mess, messid] = mkdir(output_dir);
    if s==0
        %creation failed, return with error message
        disp('Failed to create the protrusion directory');
        OK = 0;
        return;
    else
        if handles.batch_processing
            % added for batch processing, Sam, 08/22/2008
            fprintf(1,['Status: Output Directory ' output_dir ' has been created\n']);
        else
            set(handles.statusbar, 'String', ['Status: Output Directory ' output_dir ' has been created']);
        end
    end
else
    %    user_response = yesnodialog('Title','Confirm Overwrite!','String','Do you want to overwrite the results?');
    user_response = 'Yes';
    switch user_response
        case 'No'
            % take no action
            OK = 0;
            return;
        case 'Yes'
            % Prepare to close GUI application window
            [s, mess, messid] = mkdir(output_dir);
            if s==0
                %creation failed, return with error message
                disp('Failed to create the protrusion directory');
                OK = 0;
                return
            else
                if handles.batch_processing
                    % added for batch processing, Sam, 08/22/2008
                    fprintf(1,'%s\n',['Status: Output Directory ' output_dir ' has been cleaned']);
                else
                    set(handles.statusbar, 'String', ['Status: Output Directory ' output_dir ' has been cleaned']);
                end
            end
    end
end

% folder for the edge - vectors
mkdir(output_dir, 'pr_vectors');
if CC > 1
    mkdir(output_dir, 'pr_vectors_CC2');
end

QTMovieName = [output_dir filesep 'prot_movie.mov'];
MakeQTMovie('start', QTMovieName);

% create a colormap for the movie edges
cmap_edge_evolution=jet(handles.maskinfo.mask_num);
pixel_edge = cell(handles.maskinfo.mask_num, 1);
normal_matrix = cell(handles.maskinfo.mask_num-1, 1);
protrusion = cell(handles.maskinfo.mask_num-1, 1);
centroid = zeros(handles.maskinfo.mask_num, 2);
img_edge_rgb=zeros(handles.rowSize,handles.colSize,3);
cell_isolated = zeros(handles.maskinfo.mask_num, 1);

for mask_idx = 1:handles.maskinfo.mask_num
    thismask = imread([handles.directory_name filesep handles.maskinfo.maskdir(mask_idx).name]);
    L = bwlabel(thismask);
    s  = regionprops(L, 'Area','Centroid');
    Allarea=zeros(length(s),1);
    for j=1:length(s)
        Allarea(j) = s(j).Area;
    end
    [tmp1, tmp2] = max(Allarea);
    thismask = (L == tmp2);
    centroid(mask_idx,:) = s(tmp2).Centroid;
    % test if cell is isolated or touches the image border
    bi1=find(thismask(1:handles.rowSize,1)>0);
    bi2=find(thismask(1:handles.rowSize,handles.colSize)>0);
    bi3=find(thismask(1,1:handles.colSize)>0);
    bi4=find(thismask(handles.rowSize,1:handles.colSize)>0);

    if isempty(bi1) & isempty(bi2) & isempty(bi3) & isempty(bi4)
        cell_isolated(mask_idx) = 1;
        ISCLOSE = 1;
    else
        cell_isolated(mask_idx) = 0;
        ISCLOSE = 0;
    end

    uni = unique(cell_isolated(1:mask_idx));
    if length(uni) > 1
        keyboard;
        errordlg(['Some masks are open and others are close. Please go back for another segmentation'],'Bad Input','modal');
        return;
    end

    if cell_isolated(mask_idx)
        pixel_list = prOrientEdge(thismask, ORIENT_CELL, 1);
    else
        c = contourc(double(thismask), [0 0]);
        pixel_list = [c(1,2:end)' c(2,2:end)'];
    end

    dupPixel = find(sum(diff(pixel_list).^2,2) <0.0001);
    pixel_list(dupPixel,:) = [];

    if CC > 1 & (ISCLOSE == 1)
        
        
        if mask_idx >= 2
            row_diff = [centroid(mask_idx,:) - centroid(1,:)];
            center_diff = repmat(row_diff, [size(pixel_list,1) 1]);
            pixel_list = pixel_list - center_diff;
        end
    end

    if mask_idx > 1
        pixel_list_last = pixel_edge{mask_idx-1};
        if ~cell_isolated(mask_idx)
            if norm(pixel_list_last(1,:) - pixel_list(1,:)) < norm(pixel_list_last(1,:) - pixel_list(end,:))
                pixel_edge{mask_idx} = pixel_list;
            else
                pixel_edge{mask_idx} = flipud(pixel_list);
            end
        else
            pixel_edge{mask_idx} = pixel_list;
        end
    else
        pixel_edge{mask_idx} = pixel_list;
    end

    if ~cell_isolated(mask_idx)
        % close contour
        c = bwboundaries(thismask);
        if isCurveClockwise([c{1}(:,2), c{1}(:,1)])
            pixel_edge{mask_idx} = pixel_edge{mask_idx}(end:-1:1,:);
        end
    else
        % open contour
        if isCurveClockwise(pixel_edge{mask_idx})
            % make the pixel Clockwise in the image coordinate
            pixel_edge{mask_idx} = pixel_edge{mask_idx}(end:-1:1,:);
        end
    end

    %    figure, plot(pixel_edge{1}(:,1), pixel_edge{1}(:,2)); hold on; plot(pixel_edge{1}(1:10,1), pixel_edge{1}(1:10,2) ,'rx'); hold off; axis ij
    orig_pixel_list{mask_idx} = pixel_list;
    round_pixel_list = round(pixel_list);

    pixel_num = size(pixel_edge{mask_idx},1);
    for j=1:pixel_num
        if round_pixel_list(j,1) > 0 & round_pixel_list(j,1) <= handles.colSize & round_pixel_list(j,2) > 0 & round_pixel_list(j,2) <= handles.rowSize
            img_edge_rgb(round_pixel_list(j,2),round_pixel_list(j,1),1)=cmap_edge_evolution(mask_idx,1);
            img_edge_rgb(round_pixel_list(j,2),round_pixel_list(j,1),2)=cmap_edge_evolution(mask_idx,2);
            img_edge_rgb(round_pixel_list(j,2),round_pixel_list(j,1),3)=cmap_edge_evolution(mask_idx,3);
        end
    end
    h_QT = figure(100);
    imshow(img_edge_rgb, []);
    title('Boundary Evolution');
    text(5,10, ['#' num2str(mask_idx)], 'Color', 'white', 'FontSize', 12);
    MakeQTMovie addfigure
end
MakeQTMovie('finish');
close(h_QT);

if handles.batch_processing
    % added for batch processing, Sam, 08/22/2008
    fprintf(1,'%s\n',['Status: QT Movie ' QTMovieName ' has been generated.']);
else
    set(handles.statusbar, 'String', ['Status: QT Movie ' QTMovieName ' has been generated.']);
end

for time_idx = 2:handles.maskinfo.mask_num
    pixel_list = pixel_edge{time_idx};
    pixel_list_last = pixel_edge{time_idx-1};
    inputParam.time_idx = time_idx;
    inputParam.pixel_list_last = pixel_list_last;
    inputParam.pixel_list = pixel_list;
    inputParam.TOLERANCE = TOLERANCE;
    lastmask = imread([handles.directory_name filesep handles.maskinfo.maskdir(time_idx-1).name]);
    inputParam.maskTM1 = lastmask;
    inputParam.ISCLOSE = ISCLOSE;
    inputParam.dl_rate = handles.dl_rate;

    if handles.batch_processing
        % added for batch processing, Sam, 08/22/2008
        fprintf(1,'%s\n',['Status: Calculating protrusion between frame ' num2str(time_idx-1) ' and frame ' num2str(time_idx) '.']);
    else
        set(handles.statusbar, 'String', ...
            ['Status: Calculating protrusion between frame ' num2str(time_idx-1) ' and frame ' num2str(time_idx)]);
    end
    output = prSamProtrusion(inputParam);
    pixel_tm1_output = output.pixel_tm1_output;
    pixel_t_output = output.pixel_t_output;
    xy_normal = output.xy_normal;
    translate_output = output.translate_output;

    if ISCLOSE == 1
        if handles.BESTSTART == 0
            % Rotate the pixels so that the first segment at time t-1 is protruding to the first segment at time t
            % i.e.: starting pixel follows protrusion vector
            [tmp, tmpI] = min(sum((repmat(pixel_t_output(1,:), [size(pixel_list,1) 1]) - pixel_list).^2,2));
            if tmpI > 1
                pixel_edge{time_idx} = [pixel_list(tmpI:end,:); pixel_list(1:tmpI-1,:)];
            end
        end
    end

    pixel_list = pixel_edge{time_idx};
    if CC > 1 & (ISCLOSE == 1)
        row_diff = [centroid(time_idx-1,:) - centroid(1,:)];
        center_diff = repmat(row_diff, [size(pixel_tm1_output,1) 1]);
        pixel_tm1_output = pixel_tm1_output + center_diff;
        center_diff = repmat(row_diff, [size(pixel_list,1) 1]);
        pixel_list = pixel_list + center_diff;

        center_diff = repmat(row_diff, [size(pixel_t_output,1) 1]);
        pixel_t_output = pixel_t_output + center_diff;

        center_diff = repmat(row_diff, [size(output.spline_pixel_list,1) 1]);
        output.spline_pixel_list = output.spline_pixel_list + center_diff;
    end

    normal_matrix{time_idx-1} =[pixel_tm1_output xy_normal];
    protrusion{time_idx-1} = [pixel_tm1_output translate_output];
    pixel_list_last = pixel_tm1_output;

    pr_outputname = [output_dir filesep 'pr_vectors' filesep 'prot_vec_' num2str(time_idx-1)];
    if handles.batch_processing
        % added for batch processing, Sam, 08/22/2008
        fprintf(1,'%s\n',['Status: Generating protrusion figure ' pr_outputname '.tif']);
    else
        set(handles.statusbar, 'String', ...
            ['Status: Generating protrusion figure ' pr_outputname '.tif']);
    end

    clear h_prot_control;

    h_prot_control = figure('Visible','Off');
    % h_prot_control = figure(1);
    imshow(lastmask, []);
    hold on
    quiver(protrusion{time_idx-1}(:,1), protrusion{time_idx-1}(:,2), protrusion{time_idx-1}(:,3),protrusion{time_idx-1}(:,4), 0,'color','blue');
    quiver(normal_matrix{time_idx-1}(:,1), normal_matrix{time_idx-1}(:,2), normal_matrix{time_idx-1}(:,3),normal_matrix{time_idx-1}(:,4), 0, 'color','yellow');

%    plot(pixel_t_output(:,1), pixel_t_output(:,2),'k:');
    plot(pixel_t_output(:,1), pixel_t_output(:,2),'g');    
    plot(pixel_tm1_output(:,1), pixel_tm1_output(:,2),'r');
%    plot(output.spline_pixel_list(:,1), output.spline_pixel_list(:,2),'g');
    text(pixel_list(1,1), pixel_list(1,2)+5,num2str(time_idx),'BackgroundColor','green');
    plot(pixel_t_output(1,1), pixel_t_output(1,2),'go','MarkerFaceColor','green');
    text(pixel_list_last(1,1), pixel_list_last(1,2)+5,num2str(time_idx-1),'BackgroundColor','red');
    plot(pixel_tm1_output(1,1), pixel_tm1_output(1,2),'ro','MarkerFaceColor','red');

    x1 = min([pixel_list_last(:,1); pixel_list(:,1)])-5;
    x2 = max([pixel_list_last(:,1); pixel_list(:,1)])+5;
    y1 = min([pixel_list_last(:,2); pixel_list(:,2)])-5;
    y2 = max([pixel_list_last(:,2); pixel_list(:,2)])+5;
    axis([x1 x2 y1 y2]);
    axis ij
    title(handles.maskinfo.maskdir(time_idx-1).name);

    print(h_prot_control, [pr_outputname '.tif'],'-dtiff');
    print(h_prot_control, [pr_outputname '.eps'],'-depsc2');
    hgsave(h_prot_control,[pr_outputname '.fig']);
    close(h_prot_control);

    if (CC > 1) & (ISCLOSE == 1)
        row_diff = [centroid(time_idx,:) - centroid(time_idx-1,:)];
        center_diff = repmat(row_diff, [size(translate_output,1) 1]);
        translate_output = translate_output + center_diff;

        cc_pixel_list = pixel_t_output;
        center_diff = repmat(row_diff, [size(cc_pixel_list,1) 1]);
        cc_pixel_list = cc_pixel_list + center_diff;

        protrusion_CC2{time_idx-1} = [pixel_tm1_output  translate_output];

        pr_outputname_CC2 = [output_dir filesep 'pr_vectors_CC2' filesep 'prot_vec_' num2str(time_idx-1)];
        fprintf(1,'%s\n',['Status: Generating protrusion figure ' pr_outputname_CC2 '.tif']);

        clear h_prot_control_CC2;

        h_prot_control_CC2 = figure('Visible','Off');
        %        h_prot_control_CC2 = figure(2);
        imshow(lastmask, []);
        hold on
        quiver(protrusion_CC2{time_idx-1}(:,1), protrusion_CC2{time_idx-1}(:,2), protrusion_CC2{time_idx-1}(:,3),protrusion_CC2{time_idx-1}(:,4), 0,'color','blue');
        quiver(protrusion_CC2{time_idx-1}(:,1), protrusion_CC2{time_idx-1}(:,2), normal_matrix{time_idx-1}(:,3),normal_matrix{time_idx-1}(:,4), 0, 'color','yellow');

        plot(cc_pixel_list(:,1), cc_pixel_list(:,2),'g');
        plot(protrusion_CC2{time_idx-1}(:,1), protrusion_CC2{time_idx-1}(:,2),'r');

        %         text(cc_pixel_list(1,1), cc_pixel_list(1,2)+5,num2str(time_idx),'BackgroundColor','green');
        %         plot(cc_pixel_list(1,1), cc_pixel_list(1,2),'go','MarkerFaceColor','green');
        %         text(protrusion_CC2{time_idx-1}(1,1), protrusion_CC2{time_idx-1}(1,2)+5,num2str(time_idx-1),'BackgroundColor','red');
        %         plot(protrusion_CC2{time_idx-1}(1,1), protrusion_CC2{time_idx-1}(1,2),'ro','MarkerFaceColor','red');

        x1 = min([protrusion_CC2{time_idx-1}(:,1); cc_pixel_list(:,1)])-5;
        x2 = max([protrusion_CC2{time_idx-1}(:,1); cc_pixel_list(:,1)])+5;
        y1 = min([protrusion_CC2{time_idx-1}(:,2); cc_pixel_list(:,2)])-5;
        y2 = max([protrusion_CC2{time_idx-1}(:,2); cc_pixel_list(:,2)])+5;
        axis([x1 x2 y1 y2]);
        axis ij 
        title(handles.maskinfo.maskdir(time_idx-1).name);

        print(h_prot_control_CC2, [pr_outputname_CC2 '.tif'],'-dtiff');
        print(h_prot_control_CC2, [pr_outputname_CC2 '.eps'],'-depsc2');
        hgsave(h_prot_control_CC2,[pr_outputname_CC2 '.fig']);
        close(h_prot_control_CC2);
    end
end

save([output_dir filesep 'pixel_edge.mat'], 'pixel_edge','orig_pixel_list','handles');
if (CC > 1) & (ISCLOSE == 1)
    save([output_dir filesep 'protrusion.mat'],    'protrusion', 'protrusion_CC2');
else
    save([output_dir filesep 'protrusion.mat'],    'protrusion');
end
save([output_dir filesep 'normal_matrix.mat'], 'normal_matrix');
save([output_dir filesep 'centroid.mat'], 'centroid');

handles.ISCLOSE = ISCLOSE;
generatePRMap(handles);
generatePRMovie(handles);
fprintf(1,'%s\n',['Results are stored under ' output_dir]);
