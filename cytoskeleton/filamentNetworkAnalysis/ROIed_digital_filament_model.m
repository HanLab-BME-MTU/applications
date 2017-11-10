function [VIF_ROI_model, VIF_ROI_orientation_model] = ...
    ROIed_digital_filament_model(VIF_current_model,img_size, ROI)
% function to take the ROI confined filament model(with digitalization and
% orientation recalculation)

% transfer model to digital
[VIF_digital_model,VIF_orientation_model,~,~,~] ...
    = filament_model_to_digital_with_orientation(VIF_current_model);

% if the ROI is full area, just copy
if(mean2(double(ROI))==1)
    VIF_ROI_model = VIF_digital_model;
    VIF_ROI_orientation_model = VIF_orientation_model;
else
    % if there is defined ROI, check filament by filament.
    count=1;
    VIF_ROI_model = cell(1,1);
    VIF_ROI_orientation_model =  cell(1,1);
    
    for iF = 1 : length(VIF_digital_model)
        try
            x = VIF_digital_model{iF}(:,1);
            y = VIF_digital_model{iF}(:,2);
            ind_xy = sub2ind(img_size,y,x);
            roi_flag_array = ROI(ind_xy);
            [inROI_label, inROI_N] = bwlabel(roi_flag_array);
            for iR = 1 : inROI_N
                if (sum(double(inROI_label==iR))>2)
                    VIF_ROI_model{count} = [x(find(inROI_label==iR)) y(find(inROI_label==iR))];
                    VIF_ROI_orientation_model{count} = VIF_orientation_model{iF}(find(inROI_label==iR));
                    count = count +1 ;
                end
            end
        end
    end
end
