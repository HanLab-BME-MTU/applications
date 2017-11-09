function [ h_img_fig, h_scale_fig, h_width_fig ] = exploreOrientationScaleWidth( I, Rscale, Rwidth, theta )
%EXPLOREORIENTATIONSCALEWIDTH GUI to explore orientation and scale / width
%
% INPUT
% I - image to display upon which to select a pixel
%     (Could be the actual image analyzed or any image the same size as the
%     original including the Response, NMS, NLMS)
% Rscale - (optional) OrientationSpaceResponse at different scales
% Rwidth - (optional) OrientationScaleSpaceFilter response at different widths
% theta  - (optional) orientation to automatically steer to when moving the
%           pointer
%
% OUTPUT
% h_img_fig   - handle to the image figure
% h_scale_fig - handle to the orientation/scale figure
% h_width_fig - handle to the orientation/width figure
% 
% Basic Usage
% lamins.gui.exploreOrientationScaleWidth(I);
%
% Example Usage
% F = OrientationSpaceFilter.constructByRadialOrder(1/2/pi./chebpts(5,[1 4]),1,8);
% FF = OrientationScaleSpaceFilter.constructByRadialOrder(1/2/pi./chebpts(5,[1 4]),1,8);
% R = F*I;
% RR = FF*I;
% maxima = R(2).getRidgeOrientationLocalMaxima;
% lamins.gui.exploreOrientationScaleWidth(I,R,RR,maxima)
%
% Mark Kittisopikul, Feburary 2017 (earlier actually, but major revisions)
% Jaqaman Lab
% UT Southwestern

if(nargin < 2)
    % Do the filtering for the user, DEMO mode
    gcp;
    F = OrientationSpaceFilter.constructByRadialOrder(1/2/pi./chebpts(5,[1 4]),1,8);
    FF = OrientationScaleSpaceFilter.constructByRadialOrder(1/2/pi./chebpts(5,[1 4]),1,8);
    Rscale = F*I;
    Rwidth = FF*I;
    theta = Rscale(2).getRidgeOrientationLocalMaxima;
elseif(nargin < 3)
    % Just do the width filter
    F = [Rscale.filter];
    FF = OrientationScaleSpaceFilter.constructByRadialOrder([F.f_c],(F(1).f_c/F(1).b_f).^2,F(1).K);
    Rwidth = FF*I;
elseif(nargin < 4)
    % If theta is empty, do not snap to the orientation
    theta = [];
end

Rscale_orig = Rscale;
Rwidth_orig = Rwidth;

% Setup filter image overlay
% fparams.Kf = 1;
% Angular order
fparams.K = Rscale(1).filter.K;
% Filter image size
fparams.sz = 101;
% Half filter image size
fparams.hsz = floor(fparams.sz/2);
% Filter color, if empty, use impoint color
% fparams.color = [1 0 1];
fparams.color = [];
% Threshold for complete transparency
fparams.thresh = 0.1;
% Max alpha value. 1.0 = Opaque, 0.0 = Transparent
fparams.alpha = 0.7;

% filtim = real(fftshift(fft2(real(orientationSpace.kernel(1/2/pi/2,1/2/pi/2,fparams.K,0,fparams.sz)))));
% filtmask = mat2gray(abs(filtim));
% filtim = mat2gray(filtim);
% filtim(:,:,3) = 0;


h_img_fig = figure;
h_img = imshow(mat2gray(I));

hold on;
h_filt = imshow(zeros(fparams.sz));
% h_filt.AlphaData = filtmask;
hold off;
% h_img_pt = impoint;

menu = uicontextmenu(h_img_fig);
m_d = uimenu(menu,'Label',['K = ' num2str(Rscale(1).filter.K)],'Checked','on','Callback',@(~,~) selectK(Rscale(1).filter.K));
m_5 = uimenu(menu,'Label','K = 5','Checked','off','Callback',@(~,~) selectK(5));
m_3  = uimenu(menu,'Label','K = 3','Checked','off','Callback',@(~,~) selectK(3));
% m_0 = uimenu(menu,'Label','K = 0','Checked','off','Callback',@(~,~) selectK(0));
% m_neg = uimenu(menu,'Label','K = -0.5','Checked','off','Callback',@(~,~) selectK(-0.5));
m_custom = uimenu(menu,'Label','K = ?','Checked','off','Callback',@promptK);
% m_set = [m_d m_5 m_3 m_0 m_neg m_custom];
m_set = [m_d m_5 m_3 m_custom];
h_img_fig.UIContextMenu = menu;
h_img.UIContextMenu = menu;

h_scale_fig = figure;

h_width_fig = figure;

% h_quiver = [];


   

scaleWidthRange = 1/2/pi./[Rscale(1).filter.f_c Rscale(end).filter.f_c];
assertEqual(scaleWidthRange,1/2/pi./[Rwidth(1).filter.f_c Rwidth(end).filter.f_c]);

figure(h_img_fig);
h_img_pt = impoint(gca,size(I,1)/2,size(I,2)/2);
h_img_pt.setColor('m');


    % Setup axes for orientation/scale and orientation/width filters
    scaleWidthRangePoints = scaleWidthRange(1):0.1:scaleWidthRange(2);
    degreePoints = 0:179;
    
    % Initialized shared variables
    h_scale = [];
    h_scale_title = [];
    h_scale_scale_roots = [];
    h_scale_width_roots = [];
    
    h_width = [];
    h_width_title = [];
    h_width_scale_roots = [];
    h_width_width_roots = [];
    
    % Setup orientation/scale and orientation/width pointers
    h_scale_width_pts = setupScaleWidth();
    
    % Add call backs
    addNewPositionCallback(h_scale_width_pts(1),@drawQuiver);
    
    addNewPositionCallback(h_img_pt,@drawScaleWidth);
    
    % Draw the initial filter
    drawQuiver();


% while(true)
%     h_scale_width_pts = drawScaleWidth();
%     addNewPositionCallback(h_scale_width_pts(1),@drawQuiver);
%     wait(h_img_pt);
% end

    function [h_scale_width_pts] =  setupScaleWidth(~,~)
        % Obtain pixel of interest
        h_img_pt_pos = round(getPosition(h_img_pt));
        
        % Analyze width and scale at that pixel
        B_scale = OrientationScaleMatrix(real(squeeze(Rscale.getArraySpace(h_img_pt_pos(2),h_img_pt_pos(1)))),scaleWidthRange,[0 180]);
        B_width = OrientationScaleMatrix(real(squeeze(Rwidth.getArraySpace(h_img_pt_pos(2),h_img_pt_pos(1)))),scaleWidthRange,[0 180]);
        B_scale_roots = roots(diff(B_scale{:,scaleWidthRangePoints})).';
        B_width_roots = roots(diff(B_width{degreePoints,:})).';
        
        % Setup Scale Figure
        figure(h_scale_fig);
        h_scale = imagesc(B_scale);
        h_scale_title = title(sprintf('Scale Response at x = %d, y = %d, K = %d',h_img_pt_pos(1),h_img_pt_pos(2),fparams.K));
        hold on;
        h_scale_scale_roots = line( ...
            'XData',repmat(scaleWidthRangePoints(:),size(B_scale_roots,2),1), ...
            'YData',B_scale_roots(:), ...
            'Color','k', ...
            'LineStyle','none', ...
            'Marker','.');
        h_scale_width_roots = line( ...
            'XData',B_width_roots(:), ...
            'YData',repmat(degreePoints(:),size(B_width_roots,2),1), ...
            'Color','r', ...
            'LineStyle','none', ...
            'Marker','.');
        hold off;
        
        % Setup Width Figure
        figure(h_width_fig);
        h_width = imagesc(B_width);
        h_width_title = title(sprintf('Width Response at x = %d, y = %d, K = %d',h_img_pt_pos(1),h_img_pt_pos(2),fparams.K));
        hold on;
        xlabel('Width');
        h_width_scale_roots = line( ...
            'XData',repmat(scaleWidthRangePoints(:),size(B_scale_roots,2),1), ...
            'YData',B_scale_roots(:), ...
            'Color','k', ...
            'LineStyle','none', ...
            'Marker','.');
        h_width_width_roots = line( ...
            'XData',B_width_roots(:), ...
            'YData',repmat(degreePoints(:),size(B_width_roots,2),1), ...
            'Color','r', ...
            'LineStyle','none', ...
            'Marker','.');%         plot(scaleWidthRangePoints,B_scale_roots,'k.')
        
        hold off;

        % Setup pointers
        h_scale_width_pts = linkaxes_pointer([h_scale.Parent h_width.Parent]);
        h_scale_width_pts(1).setColor('m');
        h_scale_width_pts(2).setColor('m');
    end

    function [] =  drawScaleWidth(~,~)
        % Obtain pixel of interest
        h_img_pt_pos = round(getPosition(h_img_pt));

        % Analyze width and scale at that pixel
        B_scale = OrientationScaleMatrix(real(squeeze(Rscale.getArraySpace(h_img_pt_pos(2),h_img_pt_pos(1)))),scaleWidthRange,[0 180]);
        B_width = OrientationScaleMatrix(real(squeeze(Rwidth.getArraySpace(h_img_pt_pos(2),h_img_pt_pos(1)))),scaleWidthRange,[0 180]);
        B_scale_roots = roots(diff(B_scale{:,scaleWidthRangePoints})).';
        B_width_roots = roots(diff(B_width{degreePoints,:})).';
        
        % Update Orientation/Scale at pixel of interest
        h_scale_title.String = sprintf('Scale Response at x = %d, y = %d, K = %d',h_img_pt_pos(1),h_img_pt_pos(2),fparams.K);
        h_scale.CData = double(B_scale.sample(360,101));
        h_scale_scale_roots.XData = repmat(scaleWidthRangePoints(:),size(B_scale_roots,2),1);
        h_scale_scale_roots.YData = B_scale_roots(:);
        h_scale_width_roots.XData = B_width_roots(:);
        h_scale_width_roots.YData = repmat(degreePoints(:),size(B_width_roots,2),1);
        
        % Update Orientation/Width at pixel of interest
        h_width_title.String = sprintf('Width Response at x = %d, y = %d, K = %d',h_img_pt_pos(1),h_img_pt_pos(2),fparams.K);
        h_width.CData = double(B_width.sample(360,101));
        h_width_scale_roots.XData = repmat(scaleWidthRangePoints(:),size(B_scale_roots,2),1);
        h_width_scale_roots.YData = B_scale_roots(:);
        h_width_width_roots.XData = B_width_roots(:);
        h_width_width_roots.YData = repmat(degreePoints(:),size(B_width_roots,2),1);
        
        % If orientation map given, use it to steer when the image pointer
        % is moved
        if(~isempty(theta) && ~isnan(theta(h_img_pt_pos(2),h_img_pt_pos(1))))
            h_scale_width_pts_pos = getPosition(h_scale_width_pts(1));
            h_scale_width_pts(1).setPosition(h_scale_width_pts_pos(1),theta(h_img_pt_pos(2),h_img_pt_pos(1))/pi*180);
        end
        
        % Update the filter pointer
        drawQuiver();
    end

    function [] = drawQuiver(~,~)
        % This used to draw an arrow (quiver), but now it draws the filter
        % of interest
        
        % Obtain location of the image pointer
        h_img_pt_pos = round(getPosition(h_img_pt));
        % Obtain location of the orientation/scale pointer
        h_scale_width_pts_pos = getPosition(h_scale_width_pts(1));
        figure(h_img_fig);
        hold on;
        % Center filter image on point of interest
        h_filt.XData = [1 fparams.sz] -fparams.hsz -1 + h_img_pt_pos(1);
        h_filt.YData = [1 fparams.sz] -fparams.hsz -1 + h_img_pt_pos(2);
        
        % Get filter for current orientation and scale of interest
        filtim = real(fftshift(fft2(real(orientationSpace.kernel(1/2/pi/h_scale_width_pts_pos(1),1/2/pi/h_scale_width_pts_pos(1),fparams.K,h_scale_width_pts_pos(2)/180*pi,fparams.sz)))));
        filtmask = mat2gray(abs(filtim));
        filtim = mat2gray(filtim);
        if(isempty(fparams.color))
            filtim = bsxfun(@times,filtim,shiftdim(getColor(h_img_pt),-1));
        else
            filtim = bsxfun(@times,filtim,shiftdim(fparams.color,-1));
        end
%         filtim(:,:,3) = 0;
        h_filt.CData = filtim;
        h_filt.AlphaData = double(filtmask > fparams.thresh)*fparams.alpha;

%         if(~isempty(h_quiver) && isvalid(h_quiver))
%             delete(h_quiver)
%         end
%         x = h_img_pt_pos(1);
%         y = h_img_pt_pos(2);
%         u = -sind(h_scale_width_pts_pos(2));
%         v = cosd(h_scale_width_pts_pos(2));
%         scale = 10;
%         h_quiver = quiver(x-u/2*scale,y-v/2*scale,u,v,scale);
%         set(h_quiver,'LineWidth',h_scale_width_pts_pos(1));
%         hold off;
    end

    function selectK(K)
        assert(K <= Rscale_orig(1).filter.K,'K must be less than or equal to the original angular order');
        fparams.K = K;
        set(m_set,'Checked','off');
        switch(K)
            case Rscale_orig(1).filter.K
                m_d.Checked = 'on';
            case 5
                m_5.Checked = 'on';
            case 3
                m_3.Checked = 'on';
%             case 0
%                 m_0.Checked = 'on';
%             case -0.5
%                 m_neg.Checked = 'on';
            otherwise
                m_custom.Checked = 'on';
        end
                
        if(K == Rscale_orig(1).filter.K)
            Rscale = Rscale_orig;
            Rwidth = Rwidth_orig;
        else
            Rscale = Rscale_orig.getResponseAtOrderFT(K);
            Rwidth = Rwidth_orig.getResponseAtOrderFT(K);
        end
        drawScaleWidth();
    end
    function promptK(~,~)
        answer = inputdlg('Enter new K value:');
        answer = str2double(answer);
        if(~isnan(answer))
            selectK(answer);
        else
            warning('Could not parse K');
        end
    end

end

