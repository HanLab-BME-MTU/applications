
set(accept1,'Enable','on');
set(lim1,'Enable','on');
set(lim2,'Enable','on');
set(lim3,'Enable','on');
set(lim4,'Enable','on');
set(lim5,'Enable','on');
set(lim6,'Enable','on');
set(back1_name,'Enable','on');
set(back2_name,'Enable','on');


uiwait;

back_flag

if back_flag == 0;
    background_top = 0;
    f3_top=f2_top;
else
    background_top = imopen(f2_top,strel('disk',15));
    f2a_top=imsubtract(f2_top,background_top);
    f2b_top=imsubtract(imadd(f2a_top,imtophat(f2a_top,se)),imbothat(f2a_top,se));
    f3_top=imadjust(f2b_top);
end

comp1(:,:,2)=f3_top(:,:);
f3a_top=imrotate(f3_top,rot_angle);
f3b_top=imcrop(f3a_top,[p1 offset]);

imshow(f3b_top);

% set(title1,'String','New Top Post Image');

if back_flag == 0;
    background_bot = 0;
    f3_bot=f2_bot;
else
    background_bot = imopen(f2_bot,strel('disk',15));
    f2a_bot=imsubtract(f2_bot,background_bot);
    f2b_bot=imsubtract(imadd(f2a_bot,imtophat(f2a_bot,se)),imbothat(f2a_bot,se));
    f3_bot=imadjust(f2b_bot);
end

f3a_bot=imrotate(f3_bot,rot_angle);
[f3b_bot R]=imcrop(f3a_bot,[p1 offset]);

imshow(f3b_bot);

% set(title1,'String','New Bottom Post Image');

I2=f3b_top;
J2=f3b_bot;

% Find and label centers

m = round(R(2));
n = round(R(1));
HW = [R(4) R(3)];
hw = ceil(HW/ROIstep);

I4 = I2;
J4 = J2;

I7 = zeros(size(I2));
J7 = zeros(size(I2));

fprintf('\nImage I \t( M , N )\t Image J\n')
fprintf('------- \t---------\t -------\n')
i = 0;
k1 = 0;
p = 0;
j_1=1;
j_2=1;
bad_posts_top=[];
bad_posts_bot=[];


for M = 1:ROIstep:(hw(1)*ROIstep+1-ROIstep);
   i = i + 1;
   j = 0;
   for N = 1:ROIstep:(hw(2)*ROIstep+1-ROIstep);
      j = j+1;
      
      bad_post_flagA = 0;
      bad_post_flagB = 0;
      manual_flagA = 0;
      manual_flagB = 0;
      
      %%%%%%%%%%%%%%%%%%%%%
      % Local Threshold I %
      %%%%%%%%%%%%%%%%%%%%%
      
    
      A = I4(M:M+(ROIstep-1),N:N+(ROIstep-1));
      eccA=1;
      areaA = 6400; % starting area
      

     
      alp = graythresh(A); % beginning threshold value

         A_bw=im2bw(A,alp);
         A_bw2=bwmorph(A_bw,'clean');
         A_bw3=bwfill(A_bw2,'holes',8);
         A_bw4=bwmorph(A_bw3,'majority');
         A_bw5=bwmorph(A_bw4,'clean');
         regA=regionprops(bwlabel(A_bw5),'Area','Centroid','Eccentricity','Image');
         areaA=sum(cat(1,regA.Area));
         [maxA indA]=max(cat(1,regA.Area)); 
         eccA=regA(indA).Eccentricity;
         centA=(cat(1,regA.Centroid));

         
         if or(or(isempty(areaA), areaA<area_min),or(areaA>area_max,eccA>ecc_max));
            image_proc_suppA;
         end

      
     final_areaA(i,j)=areaA;
     final_eccA(i,j)=eccA;

     
      %%%%%%%%%%%%%%%%%%%%%
      % Local Threshold J %
      %%%%%%%%%%%%%%%%%%%%%
      B = J4(M:M+(ROIstep-1),N:N+(ROIstep-1));
      eccB=1;
      areaB = 6400;
      
      beta = graythresh(B); % beginning threshold value

         B_bw=im2bw(B,beta);
         B_bw2=bwmorph(B_bw,'clean');
         B_bw3=bwfill(B_bw2,'holes',8);
         B_bw4=bwmorph(B_bw3,'majority');
         B_bw5=bwmorph(B_bw4,'clean');
         regB=regionprops(bwlabel(B_bw5),'Area','Centroid','Eccentricity','Image'); 
         areaB= sum(cat(1,regB.Area));
         [maxB indB]=max(cat(1,regB.Area));
         eccB=regB(indB).Eccentricity;
         centB=(cat(1,regB.Centroid));

         
         if or(or(isempty(areaB), areaB<area_min),or(areaB>area_max,eccB>ecc_max));
            image_proc_suppB;
         end

      final_areaB(i,j)=areaB;
      final_eccB(i,j)=eccB;
    

        I6=A_bw5;
        I6a(:,:,i,j)=A_bw5;
        J6=B_bw5;
        J6a(:,:,i,j)=B_bw5;

        I7(M:M+(ROIstep-1),N:N+(ROIstep-1)) = I6;
        J7(M:M+(ROIstep-1),N:N+(ROIstep-1)) = J6;
      

          Ix(i,j) = N+centA(indA,1);
          Iy(i,j) = M+centA(indA,2);

           Jx(i,j) = N+centB(indB,1);
           Jy(i,j) = M+centB(indB,2);
         
         p = p+1;
         subplot(1,2,1),imshow(I6), title('Top Image')
         subplot(1,2,2),imshow(J6), title('Bottom Image')
         MOV(p)=getframe;

      k1 = k1 + 1;
      LongIx(k1) = Ix(i,j);
      LongIy(k1) = Iy(i,j);
      LongJx(k1) = Jx(i,j);
      LongJy(k1) = Jy(i,j);
      
      
   end
end

IIa = bwmorph(I7,'remove');
II(:,:,1)=I4;
II(:,:,2)=IIa*65535;
II(:,:,3)=0;

JJa = bwmorph(J7,'remove');
JJ(:,:,1)=J4;
JJ(:,:,2)=JJa*65535;
JJ(:,:,3)=0;


ix = LongIx;
iy = LongIy;
jx = LongJx;
jy = LongJy;

centr_top.Centroid=[transpose(ix) transpose(iy)];

centr_bot.Centroid=[transpose(jx) transpose(jy)];

%Eliminate centers near edge of image

 for i=length(centr_top):-1:1;
     if or(or(centr_top(i).Centroid(1) < 10,centr_top(i).Centroid(2) < 10),or(centr_top(i).Centroid(1)>(size(I2,2)-10),centr_top(i).Centroid(2)>(size(I2,1)-10)));
         centr_top(i)=[];
end
end
 
 for i=length(centr_bot):-1:1;
     if or(or(centr_bot(i).Centroid(1) < 10,centr_bot(i).Centroid(2) < 10),or(centr_bot(i).Centroid(1)>(size(J2,2)-10),centr_bot(i).Centroid(2)>(size(J2,1)-10)));
         centr_bot(i)=[];
end
end

centr2_top=cat(1,centr_top.Centroid);

centr2_bot=cat(1,centr_bot.Centroid);

% Calculate number of posts

nposts_top=length(centr2_top);

nposts_bot=length(centr2_bot);

y2_top=I2;
y2_bot=J2;

show_posts_top_w_NS;

uiwait(msgbox('Labeled Top Image'));

show_posts_bot_w_NS;

uiwait(msgbox('Labeled Bottom Image'));

rerun=questdlg('Rerun Data Centroid Analysis?');



 
