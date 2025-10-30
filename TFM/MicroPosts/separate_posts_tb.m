% CAL 2/5/03
% Updated 2/20/03
% Last updated 8/04/08 SH
%
% Separate posts into three categories: posts with cell to be analyzed, unoccupied posts, and posts touching other cells
% (Note: separate_posts_tb separates deflections, while separate_posts_undef_tb separates arrays of cell centers




hold off;

actin2a=imrotate(actin2,rot_angle);
actin2b=imcrop(actin2a,[p1 offset]);

comp4(:,:,2)=actin2b(:,:);

comp4(:,:,1)=f3b_top(:,:);

comp4(:,:,3)=f3b_bot(:,:);

post_actin_h1=imshow(comp4);
maximize(gcf);

hold on;

uiwait(msgbox('Click on any posts that are degenerate or are touching other cells. Press <enter> when done.','Excluded Posts','modal'));

[degen_x degen_y]=ginput;

plot(degen_x, degen_y, 'y+');

if and(isempty(degen_x),isempty(degen_y));
    degen_ok=questdlg('You have indicated that there are no collapsed posts or posts touching other cells. Is this correct?','Excluded Posts');
else
    degen_ok=questdlg('Are selected posts (yellow +) correct?','Excluded Posts');
end

while degen_ok(1:2)=='No';
    imshow(comp4);
    maximize(gcf);
    degen_x=[];
    degen_y=[];
    hold on;
    [degen_x degen_y]=ginput;
    plot(degen_x, degen_y, 'y+');
    if and(isempty(degen_x),isempty(degen_y));
        degen_ok=questdlg('You have indicated that there are no collapsed posts or posts touching other cells. Is this correct?','Excluded Posts');
    else
        degen_ok=questdlg('Are selected posts (yellow +) correct?','Excluded Posts');
    end
end

%produce only deflections of degenerated posts with other members zero
[degen_posts_x,degen_posts_y]=compare_clicks2(degen_x,degen_y,actual_x_bot,actual_y_bot,delta_x,delta_y);


% also include posts determined to be "bad" from image processing:

if isempty(bad_posts_top);
else
    for i=1:size(bad_posts_top,1);
        degen_posts_x(bad_posts_top(i,1),bad_posts_top(i,2))=delta_x(bad_posts_top(i,1),bad_posts_top(i,2));
        degen_posts_y(bad_posts_top(i,1),bad_posts_top(i,2))=delta_y(bad_posts_top(i,1),bad_posts_top(i,2));
    end
end

if isempty(bad_posts_bot);
else
    for i=1:size(bad_posts_bot,1);
        degen_posts_x(bad_posts_bot(i,1),bad_posts_bot(i,2))=delta_x(bad_posts_bot(i,1),bad_posts_bot(i,2));
        degen_posts_y(bad_posts_bot(i,1),bad_posts_bot(i,2))=delta_y(bad_posts_bot(i,1),bad_posts_bot(i,2));
    end
end

post_top_flag=0;
post_bot_flag=0;
actin_flag=1;
fn_flag=0;
vinc_flag=0;

%------------------------
% [boundaries_actin,L_actin]=bwboundaries(actin4b,4);
%     
% s_actin=regionprops(L_actin,'Area','Centroid','MajorAxisLength','MinorAxisLength','FilledArea','Perimeter'); %take binary image properties of actin
% 
% k=1;
% 
% for i=1:length(s_actin);
%     b=boundaries_actin{i};
%     if length(b)>(pi*3/convert);
%         post_actin_h(i)=plot(b(:,2),b(:,1),'y');
%         b_actin(k)=i;
%         k=k+1;
%     end
% end
% 
% roi_actin=zeros(size(y2_top,1),size(y2_top,2));
% 
% for i=1:k-1;
%     roi_actin=roi_actin+poly2mask(boundaries_actin{b_actin(i)}(:,2),boundaries_actin{b_actin(i)}(:,1),size(y2_top,1),size(y2_top,2));
% end
% 
% roi_actin_border=bwmorph(roi_actin,'remove');
% 
% %separating cell posts, exterior posts and interior posts
% l=0;
% cell_posts_position_x = [];
% cell_posts_position_y = [];
% for M = 1:ROIstep:(hw(1)*ROIstep+1-ROIstep);
%     l = l + 1;
%     j = 0;
%     for N = 1:ROIstep:(hw(2)*ROIstep+1-ROIstep);
%         j = j+1;
% 
%         if and(degen_posts_x(l,j)==0,degen_posts_y(l,j)==0);
% 
%             A_act = roi_actin(M:M+(ROIstep-1),N:N+(ROIstep-1)); %extracting ROI of actin image
%             A_post=I7(M:M+(ROIstep-1),N:N+(ROIstep-1)); %extracting ROI of posts
%             A_border = roi_actin_border(M:M+(ROIstep-1),N:N+(ROIstep-1)); %extracting ROI of actin border             
%             act_corr = A_act*double(double(A_post)); %1 if there is overlapped portion between actin and a post in ROI
% 
%             if sum(sum(act_corr)) ~= 0; % if there is overlapped portion between actin and a post in ROI
%                 cell_posts_position_x(l,j)=M;
%                 cell_posts_position_y(l,j)=N;
%             end
%         end
%     end 
% end
% %-------------------
% % set(title1,'String','Composite Actin Image with Vectors');
% 
% %assigning cell posts on purpose SH080408
% %showing current cell posts so far
% 
% hold on;
% plot(cat(1,cell_posts_position_y) + ROIstep/2,cat(1,cell_posts_position_x)+ ROIstep/2, 'c*')
% % plot(actual_x_top,actual_y_top,'r*')
% 
% hold on;
% uiwait(msgbox('Left click on any posts that considered as cell posts. Right click on posts that should not be regarded as cell posts. Press <enter> when done.','Cell Posts','modal'));
% cell_x=[];
% cell_y=[];
% [cell_x, cell_y, button]=ginput;
% for i=1: size(cell_x)
%     if button(i)==1
%         plot(cell_x(i), cell_y(i), 'go');
%     else
%         plot(cell_x(i), cell_y(i), 'ko');
%     end
% end
% 
% if and(isempty(cell_x),isempty(cell_y));
%     cell_ok=questdlg('You have indicated that there are no cell posts having no actin image touching. Is this correct?','Cell Posts');
% else
%     cell_ok=questdlg('Are selected posts (green o and black o) correct?','Cell Posts');
% end
% 
% while cell_ok(1:2)=='No';
%     imshow(comp4);
%     maximize(gcf);
%     cell_x=[];
%     cell_y=[];
%     hold on;
%     [cell_x, cell_y, button]=ginput;
%     for i=1: size(cell_x)
%         if button(i)==1
%             plot(cell_x(i), cell_y(i), 'go');
%         else
%             plot(cell_x(i), cell_y(i), 'ko');
%         end
%     end
% 
%     if and(isempty(cell_x),isempty(cell_y));
%         cell_ok=questdlg('You have indicated that there are no cell posts having no actin image touching. Is this correct?','Cell Posts');
%     else
%         cell_ok=questdlg('Are selected posts (green o and black o) correct?','Cell Posts');
%     end
% end

%make deflection of cell_x, cell_y cell_posts_x, cell_posts_x
% for M = 1:ROIstep:(hw(1)*ROIstep+1-ROIstep);
%     l = l + 1;
%     j = 0;
%     for N = 1:ROIstep:(hw(2)*ROIstep+1-ROIstep);
%         j = j+1;
%         
% [cell_posts_x,cell_posts_y]=compare_clicks2(cell_x,cell_y,actual_x_bot,actual_y_bot,delta_x,delta_y);
%-------------------------------------------------------------------------


draw_cell_tb2;%inspected by SH

for i=1:row;
    for j=1:col;
        border_posts_x(i,j)=0;
        border_posts_y(i,j)=0;
    end
end

for i=1:row;
    for j=1:col;
        if and(and(cell_posts_x(i,j)==0,degen_posts_x(i,j)==0),and(cell_posts_y(i,j)==0,degen_posts_x(i,j)==0));
            border_posts_x(i,j)=delta_x(i,j);
            border_posts_y(i,j)=delta_y(i,j);
        end
    end
end


