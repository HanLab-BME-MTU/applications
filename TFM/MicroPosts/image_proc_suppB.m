l=2;

beta = 0.6*graythresh(B); % beginning threshold value
     
       for l2=1:5;
        delta_areaB(i,j,l2)=1;
        delta_centB(i,j,l2)=10;
        delta_eccB(i,j,l2)=1;
       end
     
      betastep = graythresh(B)/100; % step
      
     while or(or(or(delta_centB(i,j,l-1) > cent_lim, delta_areaB(i,j,l-1) > area_lim),or(areaB > area_max, delta_eccB(i,j,l-1) > ecc_lim)),eccB > ecc_max); %or(eccA > eccmax, areaA > 1250);  
         beta = beta + betastep;
         B_bw=im2bw(B,beta);
         B_bw2=bwmorph(B_bw,'clean');
         B_bw3=bwfill(B_bw2,'holes',8);
         B_bw4=bwmorph(B_bw3,'majority');
         B_bw5=bwmorph(B_bw4,'clean');
         B_bw6 = bwareaopen(B_bw5,30);
         regB=regionprops(bwlabel(B_bw6),'Eccentricity','Area','Centroid','Orientation'); 
         %finding centroid for each region 073008 SH--------------
         centB_each = cat(1,regB.Centroid);
         distB_each = sqrt((ROIcent(1,1)-centB_each(:,1)).^2+(ROIcent(1,2)-centB_each(:,2)).^2);
         %comparing each distance for each indexed area 073008 SH
         [minB indB] = min(distB_each);
         areaB= regB(indB).Area;
         %--------------------------------------------
         eccB=regB(indB).Eccentricity;
         centB=regB(indB).Centroid;

%          %checking roundness
%          [B,L] = bwboundaries(areaB,'noholes');
%          stats = regionprops(L,'Area','Centroid');
%          threshold = 0.94;
%          % loop over the boundaries
%             for k = 1:length(B)
% 
%               % obtain (X,Y) boundary coordinates corresponding to label 'k'
%               boundary = B{k};
%               % compute a simple estimate of the object's perimeter
%               delta_sq = diff(boundary).^2;
%               perimeter = sum(sqrt(sum(delta_sq,2)));
%               % obtain the area calculation corresponding to label 'k'
%               area = stats(k).Area;
%               % compute the roundness metric
%               metric = 4*pi*area/perimeter^2;
% %               % display the results
% %               metric_string = sprintf('%2.2f',metric);
%               % mark objects above the threshold with a black circle
%               if metric > threshold
%                 centroid = stats(k).Centroid;
%                 plot(centroid(1),centroid(2),'ko');
%               end
%             end
         
         if or(or(isempty(areaB),beta > 0.950),eccB>0.5)
             if bad_post_flagB == 0
                 bad_post_flagB = 1;
                 beta = 0;
                 l=2;
                 indB=1;
                 eccB=1;
                 areaB=6400;
                 centB=[1 1];
             else
                image_proc_supp2B;
                break
            end
            
         elseif areaB < area_min;
             if bad_post_flagB == 0;
                 bad_post_flagB = 1;
                 beta = 0;
                 l=2;
                 indB=1;
                 eccB=1;
                 areaB=6400;
                 centB=[1 1];
             else             
                image_proc_supp2B;
                break
            end
         end
         
         all_centB(i,j,l,1)=centB(1);%centB(indB,:);
         all_centB(i,j,l,2)=centB(2);
         all_areaB(i,j,l)=areaB;
         all_eccB(i,j,l)=eccB;
         
         if l > 5;
            delta_centB(i,j,l)=sqrt((all_centB(i,j,l,1)-all_centB(i,j,l-4,1))^2+(all_centB(i,j,l,2)-all_centB(i,j,l-4,2))^2);
            delta_areaB(i,j,l)=sqrt((all_areaB(i,j,l)-all_areaB(i,j,l-4))^2)/all_areaB(i,j,l);
            delta_eccB(i,j,l)=sqrt((all_eccB(i,j,l)-all_eccB(i,j,l-4))^2)/all_eccB(i,j,l);
        end

         l=l+1;
      
      end