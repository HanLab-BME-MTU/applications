function vim_screen_network_features(labelMaskNucleus,VIF_current_seg,current_img, nms)
% function to calculate network features specifically for vim screen where
% the nucleus segmentation is present
% input:

% VIF_current_seg    the filament segmented
% current_img        the intensity(usually is image-flattened results, but
%                       in the new version, this is no per movie
%                       normalization)
% MAX_st_res         the  steerable filtering results
% nms                the nms version of ST
% 


Mask = labelMaskNucleus;

[LabelMask,NoRegion] = bwlabel(Mask);








NewLabel = randperm(NoRegion);
NewLabelMask = LabelMask;
for i = 1 :NoRegion
    NewLabelMask(LabelMask==i) = NewLabel(i);
end

[DistMask,IDX] = bwdist(Mask);
RegionMask = NewLabelMask;
RegionMask(:) = NewLabelMask(IDX);

profileIntensity
for  i = 1 :NoRegion
    region_thisCell = RegionMask==i;
    nucleus_thisCell = NewLabelMask==i;
    
    % dilating the nucleus to get the profile of the vim/mt
    % intensity/filamentdensity  as the distance increase from 0( close to
    % nucleus) to max(a quarter of the image size(1))
    
    for iDis = 1 :
    
    
    
    
end



[indy,indx] = find(Mask>0);
 h1 = figure;imagesc(RegionMask);
axis image;
 saveas(h1,'region_section_0.jpg');
hold on;
plot(indx,indy,'w.');
saveas(h1,'region_section_1.jpg');

for iL = 5 : 10:1000
    no_plot=0;
for i = 1 :NoRegion
    [indy,indx] = find(NewLabelMask==i);

    centerx = mean(indx);
    centery = mean(indy);
    
    for a = 0: 2*pi/6: 2*pi-2*pi/6
        x = centerx+iL*cos(a);
        y = centery+iL*sin(a);
        if(round(x)<=0||round(y)<=0||round(x)>size(NewLabelMask,2)||round(y)>size(NewLabelMask,1))
            continue;
        end
        if(NewLabelMask(round(y),round(x))==i)
            continue;
        end
        
        if(RegionMask(round(y),round(x))==i)
           plot(x,y,'k.'); 
        no_plot = no_plot+1;
        end
    end
    
end
saveas(h1,['region_section_',num2str(iL),'.jpg']);

if(no_plot==0)
    break;
end
end