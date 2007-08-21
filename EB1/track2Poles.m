function track2Poles

[fileName,dirName] = uigetfile('*.tif','Choose a .tif file');
I = imread([dirName,filesep,fileName]);
[xmax,ymax] = size(I);
load([dirName(1:end-8),'\point_files\config001_6p00_track_bidir.mat']);
LifeTime = 7;
indx = find( [tracks.len] >= LifeTime);
traj = tracks(indx);
leIndx = length(indx);
accum = zeros(xmax,ymax);
for i = 1:leIndx
    traj(i).endID = traj(i).startID + traj(i).len - 1;
end
poCe = 1;
s = 3;
strg=sprintf('%%.%dd',s);
for slide = 0:10:110
    b = slide;
    e = 10+slide;
    im = (b+e)/2; % get the middle image nb
    indxStr=sprintf(strg,im);
    I = imread([dirName,fileName(1:end-7),indxStr,'.tif']);
    % select tracks which fall within the 10 frame window
    indxT = find([traj.endID]>b & [traj.startID]<=e);
    traj_aux = traj(indxT);
    leIndx = length(traj_aux);
    for i = 1:leIndx
        a2(i,:)=traj_aux(i).points(end,2:-1:1); %y,x to x,y
        a1(i,:)=traj_aux(i).points(1,2:-1:1); %y,x to x,y
    end
    k = 0;
    for i = 1:leIndx
        for j = i+1:leIndx
            a = ((a2(j,1)-a1(j,1))*(a1(i,2)-a1(j,2))+(a2(j,2)-a1(j,2))*(a1(j,1)-a1(i,1)));
            b = ((a2(i,1)-a1(i,1))*(a2(j,2)-a1(j,2))-(a2(i,2)-a1(i,2))*(a2(j,1)-a1(j,1)));
            if b == 0
                b = 0.001;
            end
            alpha = a/b;
            x = round(alpha*(a2(i,1)-a1(i,1))+a1(i,1));
            y = round(alpha*(a2(i,2)-a1(i,2))+a1(i,2));
            if x > 0 && x < xmax && y > 0 && y < ymax
                k = k + 1;
                p(k,1:2) =[x y];
                accum(p(k,1),p(k,2)) = accum(p(k,1),p(k,2))+1;
            end
        end
    end
    %     figure,imshow(I,[])
    %     hold on
    %     plot(p(:,2),p(:,1),'r*')
    sigma = 5;M = 0; MM1=0; MM2 = 0;
    while (5*M)>=MM2
        accum = gauss2D(accum,sigma);
        locm=locmax2D(accum,[5 5]);
        [mx my] = find(locm);

        %         figure,imshow(accum,[])
        %         hold on
        %         plot(my,mx,'r*')
        soLocm = sort(locm(:));
        soLocm = soLocm(find(soLocm));

        %         figure,plot(soLocm(1:end))
        MM1 = soLocm(end);
        MM2 = soLocm(end-1);
        if length(soLocm)>1
            M = soLocm(end-2);
        else
            M = 0
        end
        sigma = sigma + 5 ;
    end
    figure,imshow(I,[])
    hold on
    plot(my,mx,'r*')
    [mxi(poCe) myi(poCe)]=find(accum==soLocm(end));
    plot(myi,mxi,'b-')
    plot(myi(end),mxi(end),'b*')
    poCe = poCe + 1;
end