        close all;
        clear fh;
        fh = figure;
%       manual_flagB = 0;
      
        beta = graythresh(B); % beginning threshold value

        B_bw=im2bw(B,beta);
        B_bw2=bwmorph(B_bw,'clean');
        B_bw3=bwfill(B_bw2,'holes',8);
        B_bw4=bwmorph(B_bw3,'majority');
        B_bw5=bwmorph(B_bw4,'clean');

        sp1 = subplot(2,2,1); imshow(I4); title('overall top view'); rectangle('Position',[N-1,M-1,ROIstep,ROIstep],'EdgeColor','y');
        sp2 = subplot(2,2,2); imshow(B); title(strcat(['i:' num2str(i) '   j:' num2str(j)]));
        sp3 = subplot(2,2,3); imshow(J4); title('overall bottom view(current)'); rectangle('Position',[N-1,M-1,ROIstep,ROIstep],'EdgeColor','y');
        sp4 = subplot(2,2,4); imshow(B_bw5); title('Note that this is a BOTTOM image', 'Color','r');
        maximize(gcf);

        hui_m=uicontrol(fh,'style','slider','Tag','ThreshSlide','Value',1,'min',0,'max',1/alp,'Position',[500/1152 20/800 200/1152 20/800],'Units','normalized','Callback','m04');
        %adding editing box SH120308
        hui_m05=uicontrol(fh,'style','edit','String',num2str(get(hui_m,'Value')),'Position',[700/1152 20/800 50/1152 20/800],'Units','normalized','Callback','m04_textbox');
        %-----------end--------------------
        hui_m1=uicontrol(fh,'Style', 'text','Tag','ThreshLabel','String', strcat('Threshold (0-',num2str(1/alp),'1):'),'Position',[500/1152 40/800 200/1152 20/800],'Units','normalized');
        hui_m2=uicontrol(fh,'style','pushbutton','String','OK','Position',[800/1152 20/800 40/1152 20/800],'Units','normalized','Enable','on','Callback','accept4');
        hui_m3=uicontrol(fh,'style','pushbutton','String','Manually Assign Centroid','Position',[900/1152 20/800 150/1152 40/800],'Units','normalized','Enable','on','Callback','pick_centerB');
        hui_m4=uicontrol(fh,'style','pushbutton','String','Reset','Position',[400/1152 30/800 40/1152 20/800],'Units','normalized','Enable','on','Callback','reset_postB');

        %default setting - SH120308
        number_errors = 0;
        slider.val = 1;
        guidata(fh,slider);
        %----------SH end------------------
        uiwait;

        set(hui_m3,'Visible','off');
        set(hui_m4,'Visible','off');        
        set(hui_m05,'Visible','off');  

         if manual_flagB == 0;            
            regB=regionprops(bwlabel(B_bw5),'Eccentricity','Area','Centroid','Orientation'); 
%             areaB= sum(cat(1,regB.Area));
%             [maxB indB]=max(cat(1,regB.Area));
%             eccB=regB(indB).Eccentricity;
%             centB=(cat(1,regB.Centroid));
            
            centB_each = cat(1,regB.Centroid);
            distB_each = sqrt((ROIcent(1,1)-centB_each(:,1)).^2+(ROIcent(1,2)-centB_each(:,2)).^2);
            [minB indB] = min(distB_each);
            areaB= regB(indB).Area;
            eccB=regB(indB).Eccentricity;
            centB=regB(indB).Centroid;
         
         
         else
             centB=[man_x man_y];
             eccB=1;
             areaB=0;
             B_bw5=ones(ROIstep,ROIstep);
             man_posts_bot(i,j)=1;
         end
         
         final_areaB(i,j)=areaB;
         final_eccB(i,j)=eccB;