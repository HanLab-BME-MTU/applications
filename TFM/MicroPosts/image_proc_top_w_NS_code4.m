% CAL 1/21/03
%

% Process image: read in and do background subtraction



f_top=imread(file1); % file1 = top image
% f2_top=imread(file1); % file1 = top image

% f2_top=imadjust(f_top,[double(min(min(f_top)))/2^16 double(max(max(f_top)))/2^16],[0 1]);
f2_top=imadjust(f_top);

if back_flag == 0; %when is it zero?
    background_top = 0;
    f3_top=f2_top;
else
    background_top = imopen(f2_top,strel('disk',15));
    f2a_top=imsubtract(f2_top,background_top);
%    f2b_top=imsubtract(imadd(f2a_top,imtophat(f2a_top,se)),imbothat(f2a_top,se));
%    f3_top=imadjust(f2b_top);
    f3_top=f2a_top;
end


% build composite of actin and top images

comp1(:,:,1)=f3_top(:,:);

comp1(:,:,2)=actin2(:,:);

comp1(:,:,3)=0;

% Determine cropping and rotation


subplot(1,1,1);%back into a single screeen

imshow(comp1);
maximize(gcf);

uiwait(msgbox('Click on upper left post to be analyzed and upper right post to be analyzed to set rotation angle. Click OK to continue.','mage Rotation','modal'));

hold on;

[rot_y1, rot_x1]=ginput(1);

plot(rot_y1,rot_x1,'w+');

[rot_y2, rot_x2]=ginput(1);

plot(rot_y2,rot_x2,'w+');

rot_ok=questdlg('Are rotation points (blue) correct?','Rotation Points');

while rot_ok(1:2)=='No';
    imshow(comp1);
    maximize(gcf);
    hold on;
    [rot_y1 rot_x1]=ginput(1);
    plot(rot_y1,rot_x1,'w+');
    [rot_y2 rot_x2]=ginput(1);
    plot(rot_y2,rot_x2,'w+');
    rot_ok=questdlg('Are rotation points (blue) correct?','Rotation Points');
end

posts_between=inputdlg(['How many posts between blue rotation points?'; '(including posts marked by rotation points):';'                                            '],'Post Spacing');

y_spacing=sqrt((rot_x2-rot_x1)^2+(rot_y2-rot_y1)^2)/(str2num(posts_between{1})-1);

ROIstep=floor(y_spacing);

convert=pin_spacing/y_spacing; 

x_spacing=y_spacing;


rot_angle=-atan2((rot_x1-rot_x2),(rot_y2-rot_y1))*180/pi;

f3a_top=imrotate(f3_top,rot_angle);

actin4a=imrotate(actin4,rot_angle);

comp1a=imrotate(comp1,rot_angle);

hold off;

imshow(comp1a);
maximize(gcf);

uiwait(msgbox(['Drag Rectangle to Crop Image. Set the Crop Box outside of the posts to be analyzed such that it is not touching any posts. Crop box will automatically resize to units of ' num2str(ROIstep) ' pixels by ' num2str(ROIstep) ' pixels. Click OK to Continue.'],'Image Cropping','modal'));
%need to upgrade this for centers to be the center of ROI 080208 SH
    k = waitforbuttonpress;
    point1 = get(gca,'CurrentPoint');    % button down detected
    finalRect = rbbox;                   % return figure units
    point2 = get(gca,'CurrentPoint');    % button up detected
    point1 = point1(1,1:2);              % extract x and y
    point2 = point2(1,1:2);
    p1 = min(point1,point2);             % calculate locations
    o1 = abs(point1-point2);         % and dimensions
    offset=ceil(o1/ROIstep)*ROIstep;
    x = [p1(1) p1(1)+offset(1) p1(1)+offset(1) p1(1) p1(1)];
    y = [p1(2) p1(2) p1(2)+offset(2) p1(2)+offset(2) p1(2)];
    hold on
    axis manual
    plot(x,y,'y','linewidth',5)
    
 crop_ok=questdlg('Is Crop Area (yellow) ok?','Crop Area');
 
while crop_ok(1:2)=='No'; 
    imshow(comp1a);
    maximize(gcf);
    k = waitforbuttonpress;
    point1 = get(gca,'CurrentPoint');    % button down detected
    finalRect = rbbox;                   % return figure units
    point2 = get(gca,'CurrentPoint');    % button up detected
    point1 = point1(1,1:2);              % extract x and y
    point2 = point2(1,1:2);
    p1 = min(point1,point2);             % calculate locations
    o1 = abs(point1-point2);         % and dimensions
    offset=ceil(o1/ROIstep)*ROIstep;
    x = [p1(1) p1(1)+offset(1) p1(1)+offset(1) p1(1) p1(1)];
    y = [p1(2) p1(2) p1(2)+offset(2) p1(2)+offset(2) p1(2)];
    hold on
    axis manual
    plot(x,y,'y','linewidth',5)
    crop_ok=questdlg('Is Crop Area (yellow) ok?','Crop Area');
end
    
f3b_top=imcrop(f3a_top,[p1 offset]);
actin4b=imcrop(actin4a,[p1 offset]);
comp1b=imcrop(comp1a,[p1 offset]);


%Process bottom image to get the centers of pillars

f_bot=imread(file2);
% f2_bot=imread(file2);

f2_bot=imadjust(f_bot,[double(min(min(f_bot)))/2^16 double(max(max(f_bot)))/2^16],[0 1]);
% f2_bot=imadjust(f_bot);

se=strel('disk',3);


if back_flag == 0;
    background_bot = 0;
    f3_bot=f2_bot;
else
    background_bot = imopen(f2_bot,strel('disk',15));
    f2a_bot=imsubtract(f2_bot,background_bot);
%    f2b_bot=imsubtract(imadd(f2a_bot,imtophat(f2a_bot,se)),imbothat(f2a_bot,se));
%    f3_bot=imadjust(f2b_bot);
    f3_bot=f2a_bot;
end

f3a_bot=imrotate(f3_bot,rot_angle);
[f3b_bot R]=imcrop(f3a_bot,[p1 offset]);    

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


for M = 1:ROIstep:(hw(1)*ROIstep+1-ROIstep); %row
   i = i + 1;
   j = 0;
   for N = 1:ROIstep:(hw(2)*ROIstep+1-ROIstep); %column
      j = j+1;
      
      bad_post_flagA = 0;
      bad_post_flagB = 0;
      manual_flagA = 0;
      manual_flagB = 0;
      
      %%%%%%%%%%%%%%%%%%%%%
      % Local Threshold I %
      %%%%%%%%%%%%%%%%%%%%%
      
    
      A = I4(M:M+(ROIstep-1),N:N+(ROIstep-1)); %M:rows(y-axis), N:columns(x-axis)
      eccA=1;
      areaA = 6400; % starting area(why??)
      

     
      alp = graythresh(A); % beginning threshold value

         A_bw=im2bw(A,alp);
         A_bw2=bwmorph(A_bw,'clean');
         A_bw3=bwfill(A_bw2,'holes',8);
         A_bw4=bwmorph(A_bw3,'majority');
         A_bw5=bwmorph(A_bw4,'clean');
         regA=regionprops(bwlabel(A_bw5),'Area','Centroid','Eccentricity','Image');
        
         %finding centroid for each region 073008 SH
         ROIcent = [ROIstep/2,ROIstep/2];
         centA_each = cat(1,regA.Centroid);
         distA_each = sqrt((ROIcent(1,1)-centA_each(:,1)).^2+(ROIcent(1,2)-centA_each(:,2)).^2);
         %comparing each distance for each indexed area
         [minA indA] = min(distA_each);
         %[maxA indA]=max(cat(1,regA.Area)); 
         
         areaA=regA(indA).Area;
         eccA=regA(indA).Eccentricity;
         centA=regA(indA).Centroid;

         
         if or(or(isempty(areaA), areaA<area_min),or(areaA>area_max,eccA>ecc_max));
            image_proc_suppA;
         end

      
     final_areaA(i,j)=areaA;
     final_eccA(i,j)=eccA;

     
      %%%%%%%%%%%%%%%%%%%%%
      % Local Threshold J % For bottom (SH)
      %%%%%%%%%%%%%%%%%%%%%
      B = J4(M:M+(ROIstep-1),N:N+(ROIstep-1));
      eccB=1;
      areaB = 6400;
      Bfiltered = filterGauss2D(B,2);
      if isa(B,'uint16')
          Bfiltered=uint16(Bfiltered);
      end
      beta = graythresh(Bfiltered); % beginning threshold value

         B_bw=im2bw(Bfiltered,beta);
         B_bw2=bwmorph(B_bw,'clean');
         B_bw3=bwfill(B_bw2,'holes',8);
         B_bw4=bwmorph(B_bw3,'majority');
         B_bw5=bwmorph(B_bw4,'clean');
         regB=regionprops(B_bw5,'Area','Centroid','Eccentricity','Image'); 
         %finding centroid for each region 073008 SH
         centB_each = cat(1,regB.Centroid);
         distB_each = sqrt((ROIcent(1,1)-centB_each(:,1)).^2+(ROIcent(1,2)-centB_each(:,2)).^2);
         %comparing each distance for each indexed area 073008 SH
         [minB, indB] = min(distB_each);
         areaB= regB(indB).Area;
         %[maxB indB]=max(cat(1,regB.Area));
         eccB=regB(indB).Eccentricity;
         centB=regB(indB).Centroid;

         
         if or(or(isempty(areaB), areaB<area_min),or(areaB>area_max,eccB>ecc_max))
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
      

        %saving centers into Ix,Iy (top) and Jx, Jy (bottom) SH
          Ix(i,j) = N+centA(1,1);
          Iy(i,j) = M+centA(1,2);

           Jx(i,j) = N+centB(1,1);
           Jy(i,j) = M+centB(1,2);
         
         p = p+1;
         subplot(1,2,1),imshow(I6), title('Top Image')
         hold on;
         plot(Ix(i,j),Iy(i,j),'b*');

         subplot(1,2,2),imshow(J6), title('Bottom Image')
         hold on;
         plot(Jx(i,j),Jy(i,j),'r*');
         MOV(p)=getframe;
         hold off;

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

%  for i=length(centr_top):-1:1;
%      if or(or(centr_top(i).Centroid(1) < 10,centr_top(i).Centroid(2) < 10),or(centr_top(i).Centroid(1)>(size(I2,2)-10),centr_top(i).Centroid(2)>(size(I2,1)-10)));
%          centr_top(i)=[];
% end
% end
%  
%  for i=length(centr_bot):-1:1;
%      if or(or(centr_bot(i).Centroid(1) < 10,centr_bot(i).Centroid(2) < 10),or(centr_bot(i).Centroid(1)>(size(J2,2)-10),centr_bot(i).Centroid(2)>(size(J2,1)-10)));
%          centr_bot(i)=[];
% end
% end

%final centers through lots of conversions... centr2_top, centr2_bot
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

while rerun(1:2)=='Ye'
    image_proc_top_w_NS_code4_repeat2;
end



%rearrange grids to fit form needed later in program:

y_top=double(centr2_top(:,1));
x_top=double(centr2_top(:,2));

l=length(centr2_bot);

y_bot=double(centr2_bot(:,1));
x_bot=double(centr2_bot(:,2));
 
