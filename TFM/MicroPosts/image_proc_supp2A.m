      close all;
      clear fh;
      fh = figure;
 
 
      alp=graythresh(A);
         A_bw=im2bw(A,alp);
         A_bw2=bwmorph(A_bw,'clean');
         A_bw3=bwfill(A_bw2,'holes',8);
         A_bw4=bwmorph(A_bw3,'majority');
         A_bw5=bwmorph(A_bw4,'clean');
      
         %Note: For graphic object, the first entry represents x axis
         %direction whereas first entry of image data means n-th of rows
         %which means y-axis - written by SH 120408
      subplot(2,2,1); imshow(I4); title('overall top view(current)'); rectangle('Position',[N-1,M-1,ROIstep,ROIstep],'EdgeColor','y');
      subplot(2,2,2); imshow(A); title(strcat(['i:' num2str(i) '   j:' num2str(j)]));
      subplot(2,2,3); imshow(J4); title('overall bottom view'); rectangle('Position',[N-1,M-1,ROIstep,ROIstep],'EdgeColor','y');
      subplot(2,2,4); imshow(A_bw5); title('Note that this is a TOP image','Color','r');
      maximize(gcf);
        
         hui_m=uicontrol(fh,'style','slider','Tag','ThreshSlide','Value',1,'min',0,'max',1/alp,'Position',[500/1152 20/800 200/1152 20/800],'Units','normalized','Callback','m03');
         %adding editing box SH120308
         hui_m05=uicontrol(fh,'style','edit','String',num2str(get(hui_m,'Value')),'Position',[700/1152 20/800 50/1152 20/800],'Units','normalized','Callback','m03_textbox');
         %-----------end--------------------
         hui_m1=uicontrol(fh,'Style', 'text','Tag','ThreshLabel','String', strcat('Threshold (0-',num2str(1/alp),'1):'),'Position',[500/1152 40/800 200/1152 20/800],'Units','normalized');
         hui_m2=uicontrol(fh,'style','pushbutton','String','OK','Position',[800/1152 20/800 40/1152 20/800],'Units','normalized','Enable','on','Callback','accept3');
         hui_m3=uicontrol(fh,'style','pushbutton','String','Manually Assign Centroid','Position',[900/1152 20/800 150/1152 40/800],'Units','normalized','Enable','on','Callback','pick_centerA');
         hui_m4=uicontrol(fh,'style','pushbutton','String','Reset','Position',[400/1152 30/800 40/1152 20/800],'Units','normalized','Enable','on','Callback','reset_post');

         %default setting - SH120308
         number_errors = 0;
         slider.val = 1;
         guidata(fh,slider);
         %----------SH end------------------
         
         uiwait;
         
         set(hui_m05, 'Visible','off');
         set(hui_m3,'Visible','off');
         set(hui_m4,'Visible','off');
         
         if manual_flagA == 0;
            regA=regionprops(bwlabel(A_bw5),'Eccentricity','Area','Centroid','Orientation');
            areaA=sum(cat(1,regA.Area));
            [maxA indA]=max(cat(1,regA.Area));   
            eccA=regA(indA).Eccentricity;
            centA=(cat(1,regA.Centroid));
         else
             centA=[man_y man_x];
             eccA=1;
             areaA=0;
             man_posts_top(i,j)=1;
             A_bw5=ones(ROIstep,ROIstep);

         end
         
         final_areaA(i,j)=areaA;
         final_eccA(i,j)=eccA;
 
  

     
     