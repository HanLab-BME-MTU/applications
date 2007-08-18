function trackPoles

I = imread(['C:\amatov\data\Monastrol\images\mono3001.tif']);
[xmax,ymax] = size(I);
load(['C:\amatov\data\Monastrol\point_files\config001_6p00_track_bidir.mat']);

LifeTime = 7;
indx = find( [tracks.len] >= LifeTime);
traj = tracks(indx);
leIndx = length(indx);
accum = zeros(xmax,ymax);

for i = 1:leIndx
    a2(i,:)=traj(i).points(end,2:-1:1); %y,x to x,y
    a1(i,:)=traj(i).points(1,2:-1:1); %y,x to x,y
end
k = 0;
for i = 1:leIndx
    for j = i+1:leIndx
        alpha = ((a2(j,1)-a1(j,1))*(a1(i,2)-a1(j,2))+(a2(j,2)-a1(j,2))*(a1(j,1)-a1(i,1)))/...
            ((a2(i,1)-a1(i,1))*(a2(j,2)-a1(j,2))-(a2(i,2)-a1(i,2))*(a2(j,1)-a1(j,1)));

        x = round(alpha*(a2(i,1)-a1(i,1))+a1(i,1));
        y = round(alpha*(a2(i,2)-a1(i,2))+a1(i,2));
        if x > 0 && x < xmax && y > 0 && y < ymax
            k = k + 1;
            p(k,1:2) =[x y];
            accum(p(k,1),p(k,2)) = accum(p(k,1),p(k,2))+1;
        end
    end
end

figure,imshow(accum,[])


figure,imshow(I,[])
hold on
plot(p(:,2),p(:,1),'r*')
hold off
p


% figure,imshow(ac,[])
while M<MM
  ac = gauss2D(ac,4);
% figure,imshow(ac2,[])
% ac3 = gauss2d(ac2,4);
% figure,imshow(ac3,[])
% 
% locm=locmax2D(ac3,[5 5]);
% [mx my] = find(locm);
% hold on
% plot(mx,my,'r*')
% 
% figure,imshow(ac3,[])
% hold on
% plot(my,mx,'r*')
% 
% soLocm = sort(locm(:));
% soLocm = soLocm(find(soLocm));
%figure,plot(soLocm(end-50:end))
% figure,plot(soLocm)
% size(mx)
% ans =
%    228     1
% soLocm(end)
% ans =
%     2.2762
% soLocm(end-1)
% ans =
%     0.2760
% 
% figure,imshow(I,[])
% hold on
% plot(my,mx,'r*')
% 
% [mxi myi]=find(ac3==soLocm(end))
% mxi =
%    226
% myi =
%    351
% hold on
% plot(myi,mxi,'b*')


