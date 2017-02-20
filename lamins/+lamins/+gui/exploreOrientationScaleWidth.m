function [ h ] = exploreOrientationScaleWidth( I, Rscale, Rwidth )
%EXPLOREORIENTATIONSCALEWIDTH Summary of this function goes here
%   Detailed explanation goes here

h_img_fig = figure;
h_img = imshow(I,[]);
h_img_pt = impoint;

h_scale_fig = figure;

h_width_fig = figure;

h_quiver = [];

scaleWidthRange = 1/2/pi./[Rscale(1).filter.f_c Rscale(end).filter.f_c];
assertEqual(scaleWidthRange,1/2/pi./[Rwidth(1).filter.f_c Rwidth(end).filter.f_c]);

while(true)
    h_scale_width_pts = drawScaleWidth();
    addNewPositionCallback(h_scale_width_pts(1),@drawQuiver);
    wait(h_img_pt);
end

    function [h_scale_width_pts] =  drawScaleWidth(~,~)
        h_img_pt_pos = round(getPosition(h_img_pt));
        
        B_scale = OrientationScaleMatrix(real(squeeze(Rscale.getArraySpace(h_img_pt_pos(2),h_img_pt_pos(1)))),scaleWidthRange,[0 180]);
        B_width = OrientationScaleMatrix(real(squeeze(Rwidth.getArraySpace(h_img_pt_pos(2),h_img_pt_pos(1)))),scaleWidthRange,[0 180]);
        B_scale_roots = roots(diff(B_scale{:,scaleWidthRange(1):0.1:scaleWidthRange(2)}));
        B_width_roots = roots(diff(B_width{0:179,:}));
        
        figure(h_scale_fig);
        h_scale = imagesc(B_scale);
        hold on;
        plot(scaleWidthRange(1):0.1:scaleWidthRange(2),B_scale_roots,'k.')
        plot(B_width_roots,0:179,'r.')
        hold off;
        
        figure(h_width_fig);
        h_width = imagesc(B_width);
        hold on;
        xlabel('Width');
        plot(scaleWidthRange(1):0.1:scaleWidthRange(2),B_scale_roots,'k.')
        plot(B_width_roots,0:179,'r.')
        
        hold off;

        h_scale_width_pts = linkaxes_pointer([h_scale.Parent h_width.Parent]);

    end

    function [] = drawQuiver(~,~)
        h_img_pt_pos = round(getPosition(h_img_pt));
        h_scale_width_pts_pos = getPosition(h_scale_width_pts(1));
        figure(h_img_fig);
        hold on;
        if(~isempty(h_quiver) && isvalid(h_quiver))
            delete(h_quiver)
        end
        h_quiver = quiver(h_img_pt_pos(1),h_img_pt_pos(2),-sind(h_scale_width_pts_pos(2)),cosd(h_scale_width_pts_pos(2)),10);
        set(h_quiver,'LineWidth',h_scale_width_pts_pos(1));
        hold off;
    end

end

