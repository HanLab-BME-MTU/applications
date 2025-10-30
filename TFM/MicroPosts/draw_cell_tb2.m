% Determine outline of cell and use to identify which posts cells are on








[boundaries_actin,L_actin]=bwboundaries(actin4b,4);
    
s_actin=regionprops(L_actin,'Area','Centroid','MajorAxisLength','MinorAxisLength','FilledArea','Perimeter'); %take binary image properties of actin

k=1;

for i=1:length(s_actin);
    b=boundaries_actin{i};
    if length(b)>(pi*3/convert);
        post_actin_h(i)=plot(b(:,2),b(:,1),'y');
        b_actin(k)=i;
        k=k+1;
    end
end

cell_num=0;
exterior_num=0;
interior_num=0;

%initializing
for i=1:row;
    for j=1:col;
        cell_posts_x(i,j)=0;
        cell_posts_y(i,j)=0;
        exterior_posts_x(i,j)=0;
        exterior_posts_y(i,j)=0;
        interior_posts_x(i,j)=0;
        interior_posts_y(i,j)=0;
    end
end

roi_actin=zeros(size(y2_top,1),size(y2_top,2));

for i=1:k-1;
    roi_actin=roi_actin+poly2mask(boundaries_actin{b_actin(i)}(:,2),boundaries_actin{b_actin(i)}(:,1),size(y2_top,1),size(y2_top,2));
end

roi_actin_border=bwmorph(roi_actin,'remove');

%separating cell posts, exterior posts and interior posts
l=0;
for M = 1:ROIstep:(hw(1)*ROIstep+1-ROIstep);
    l = l + 1;
    j = 0;
    for N = 1:ROIstep:(hw(2)*ROIstep+1-ROIstep);
        j = j+1;

        if and(degen_posts_x(l,j)==0,degen_posts_y(l,j)==0);

            A_act = roi_actin(M:M+(ROIstep-1),N:N+(ROIstep-1)); %extracting ROI of actin image
            A_post=I7(M:M+(ROIstep-1),N:N+(ROIstep-1)); %extracting ROI of posts
            A_border = roi_actin_border(M:M+(ROIstep-1),N:N+(ROIstep-1)); %extracting ROI of actin border             
            act_corr = A_act*double(double(A_post)); %1 if there is overlapped portion between actin and a post in ROI

            final_areaA(i,j)=0;
            final_areaB(i,j)=0;
            final_eccA(i,j)=0;
            final_eccB(i,j)=0;

            if sum(sum(act_corr)) ~= 0; % if there is overlapped portion between actin and a post in ROI
                cell_posts_x(l,j)=delta_x(l,j); %saving deflection in cell posts arrays
                cell_posts_y(l,j)=delta_y(l,j);
                cell_num=cell_num+1;%see this for inspection
                cell_posts_position_x(l,j)=M;
                cell_posts_position_y(l,j)=N;
                %else
                %cell_posts_x(j,l)=0;
                %cell_posts_y(j,l)=0;
                if sum(sum(A_border)) ~= 0;
                    exterior_posts_x(l,j)=delta_x(l,j);
                    exterior_posts_y(l,j)=delta_y(l,j);
                    exterior_num=exterior_num+1;
                else
                    interior_posts_x(l,j)=delta_x(l,j);
                    interior_posts_y(l,j)=delta_y(l,j);
                    interior_num=interior_num+1;
                end

            end
        end
    end 
end






