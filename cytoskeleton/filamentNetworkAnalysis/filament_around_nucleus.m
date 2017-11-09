function network_feature = load_MD_network_for_analysis(MD,ROI,radius)
% function to do network analysis with input MD

% input:    MD:    the loaded movieData object.
%           ROI:   the user input ROI, if not defined [], will use the whole
%                  image area
% output:   network_feature, a cell structure for each channel, each frame.
%           Each struct with field of:
%           'straightness_per_filament_pool','orientation_pixel_pool_display',...
%           'length_per_filament_pool'. 
            
% Liya Ding 2013

movie_Dir = MD.outputDirectory_;

% find all the index of diffe
avg_int_cell = cell(1,1);
sum_int_cell = cell(1,1);
matfilenames = {'M:\Work\Alexis_cell\Sept_2013\movieData.mat'};

colorarray = rand(3,20);

for iMovie = 1:1
    load(matfilenames{iMovie});
    
    avg_int_array = nan(MD.nFrames_,20);
    sum_int_array = nan(MD.nFrames_,20);
    
    for iFrame = 1:MD.nFrames_
        im = [];
        mask=[];
        im_org = double(MD.channels_(1).loadImage(iFrame));
        
        mask_org = MD.processes_{2}.loadChannelOutput(1,iFrame);

        im = imresize(im_org, [3*size(im_org,1) 3*size(im_org,2)],'bilinear');
        mask = imresize(double(mask_org), [3*size(mask_org,1) 3*size(mask_org,2)],'bilinear')>0.5;
        
        
        big_mask = imdilate(mask,ones(50,50));
        
        im = im - mean(im(big_mask==0));
        im(im<0)=0;
       
        max_int = double(mean((im(mask>0))));
        
        mask_array = cell(1,1000);
        mask_array{1}=mask;
        figure(2); imagescc(im);
        
     
        labelMask = bwlabel(mask);
        ob_prop = regionprops(labelMask,'Area','MajorAxisLength','BoundingBox');
        
        obAreas = [ob_prop.Area];
        obLongaxis = [ob_prop.MajorAxisLength];
        obBox =  [ob_prop.BoundingBox];
        box_size = min(obBox(3:4));
        step = floor(box_size/42);
        
%         avg_int_array(iFrame,j) = mean(im((mask_array{j}-mask_array{j+1})>0))/max_int;
%         sum_int_array(iFrame,j) = sum(im((mask_array{j}-mask_array{j+1})>0))/double(sum(sum(im(mask>0))));
        temp_all_masks = cell(1,1000);
        temp=mask;
        for m = 1 : 1000
            temp = imerode(temp,[0 1 0; 1 1 1; 0 1 0]);
            temp_all_masks{m} = temp;
            if(sum(sum( temp))<4)
                break;
            end
        end
        temp_all_masks{m+1} = temp>100;
        
        size_number = m;
        
        step_array = round(1:double(m)/20:m+1);
        
        for j = 1 :20
            mask_array{j+1} = temp_all_masks{step_array(j+1)};
            [indy,indx] = find((mask_array{j}-mask_array{j+1})>0);
            hold on; plot(indx, indy,'.','color',colorarray(:,j));
            avg_int_array(iFrame,j) = mean(im((mask_array{j}-mask_array{j+1})>0))/max_int;
            sum_int_array(iFrame,j) = sum(im((mask_array{j}-mask_array{j+1})>0))/double(sum(sum(im(mask>0))));            
        end
        
        avg_int_cell{iMovie} = avg_int_array;
        sum_int_cell{iMovie} = sum_int_array;
    end
end
    
figure;bar(sum_int_array');
legend({'1', '2'});
 
ylabel('Total intensity in this region in the total intenisty');
xlabel('Distance from cell boundary, in % of cell radius');
set(gca,'XTick',1:20);
set(gca,'XTickLabel',{'5','10','15','20',...
    '25','30','35','40',...
    '45','50','55','60',...
    '65','70','75','80',...
    '85','90','95','100'});
  
title('Sum intensity profiles');
axis auto;

figure;bar(avg_int_array');
legend({'1', '2'});
axis auto;
 
ylabel('Mean intensity of subregion/cell mean intensity');
xlabel('Distance from cell boundary, in % of cell radius');
set(gca,'XTick',1:20);
set(gca,'XTickLabel',{'5','10','15','20',...
    '25','30','35','40',...
    '45','50','55','60',...
    '65','70','75','80',...
    '85','90','95','100'});
title('Mean intenisty profiles');
axis auto;
