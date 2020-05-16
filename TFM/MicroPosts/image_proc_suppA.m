      k=2; % indices for looping through thresholds
          
      alp = 0.6*graythresh(A); % beginning threshold value
      
      for k2=1:5;               % dummy values needed for averaging over 5 values
        delta_areaA(i,j,k2)=1;
        delta_centA(i,j,k2)=10;
        delta_eccA(i,j,k2)=1;
      end
      
      alpstep = graythresh(A)/100; % step size
      
      while or(or(or(delta_centA(i,j,k-1) > cent_lim, delta_areaA(i,j,k-1) > area_lim),or(areaA > area_max, delta_eccA(i,j,k-1) > ecc_lim)),eccA > ecc_max); %or(eccA > eccmax, areaA > 1250);  
         alp = alp + alpstep;
         A_bw=im2bw(A,alp);
         A_bw2=bwmorph(A_bw,'clean');
         A_bw3=bwfill(A_bw2,'holes',8);
         A_bw4=bwmorph(A_bw3,'majority');
         A_bw5=bwmorph(A_bw4,'clean');
         regA=regionprops(bwlabel(A_bw5),'Eccentricity','Area','Centroid','Orientation');
         %inserted by SH-----------------------------------------
         centA_each = cat(1,regA.Centroid);
         distA_each = sqrt((ROIcent(1,1)-centA_each(:,1)).^2+(ROIcent(1,2)-centA_each(:,2)).^2);
         [minA indA] = min(distA_each);
         %-------------------------------------------------------
         %[maxA indA]=max(cat(1,regA.Area));   
         %areaA=sum(cat(1,regA(indA).Area));
         areaA=regA(indA).Area;
         eccA=regA(indA).Eccentricity;
         centA=regA(indA).Centroid;

         if or(isempty(areaA), alp >0.950)
             if bad_post_flagA == 0;
                bad_post_flagA = 1;
                alp = 0;
                k=2;
                indA=1;
                eccA=1;
                areaA=6400;
                centA=[1 1];
            else
                image_proc_supp2A;
                break
            end
            
         elseif areaA < area_min
             if bad_post_flagA == 0;
                 bad_post_flagA = 1;
                 alp=0;
                 k=2;
                 indA=1;
                 eccA=1;
                 areaA=6400;
                 centA=[1 1];
             else
                image_proc_supp2A;
                break
            end
         end
        
        
         all_centA(i,j,k,1)=centA(1);%,centA(indA,:);
         all_centA(i,j,k,2)=centA(2);
         all_areaA(i,j,k)=areaA;
         all_eccA(i,j,k)=eccA;

         if k > 5;
                delta_centA(i,j,k)=sqrt((all_centA(i,j,k,1)-all_centA(i,j,k-4,1))^2+(all_centA(i,j,k,2)-all_centA(i,j,k-4,2))^2);
                delta_areaA(i,j,k)=sqrt((all_areaA(i,j,k)-all_areaA(i,j,k-4))^2)/all_areaA(i,j,k);
                delta_eccA(i,j,k)=sqrt((all_eccA(i,j,k)-all_eccA(i,j,k-4))^2)/all_eccA(i,j,k);
         end
         k=k+1;
  
     end