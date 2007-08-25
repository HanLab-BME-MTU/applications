function track3Poles

[fileName,dirName] = uigetfile('*.tif','Choose a .tif file');
I = imread([dirName,filesep,fileName]);
[xmax,ymax] = size(I);
load([dirName(1:end-8),'\point_files\config001_5p00_track_bidir.mat']);
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
window = 2;
for begin = (LifeTime+1):window:(125-LifeTime)
    b = begin;
    e = window+begin;
    im = b-LifeTime;%round((b+e)/2); % get the middle image nb
    indxStr=sprintf(strg,im);
    I = imread([dirName,fileName(1:end-7),indxStr,'.tif']);
    % select tracks which fall within the 10 frame window
    indxT = find([traj.endID]>=b & [traj.startID]<e);
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
    sigma = 4;M = 0; MM1=0; MM2 = 0;

    accum = gauss2D(accum,sigma);
    locm=locmax2D(accum,[5 5]);
    [mx my] = find(locm);

    %                 figure,imshow(accum,[])
    %                 hold on
    %                 plot(my,mx,'r*')
    soLocm = sort(locm(:));
    soLocm = soLocm(find(soLocm));

    %                 figure,plot(soLocm(1:end))
    MM1 = soLocm(end);

    MM2 = soLocm(end-1);


    figure,imshow(I,[])
    hold on
    %     plot(my,mx,'r*')
    [x1 y1]=find(accum==soLocm(end));
    if length(soLocm) > 1
        [x2 y2]=find(accum==soLocm(end-1));
        if y1 > y2
            mxi1(poCe) = x1;
            myi1(poCe) = y1;
            mxi2(poCe) = x2;
            myi2(poCe) = y2;
        else
            mxi1(poCe) = x2;
            myi1(poCe) = y2;
            mxi2(poCe) = x1;
            myi2(poCe) = y1;
        end
    else
        mxi2(poCe) = x1 + 10;
        myi2(poCe) = y1 + 10;
        mxi1(poCe) = x1;
        myi1(poCe) = y1;
    end
    plot(myi1,mxi1,'b-')
    plot(myi1(end),mxi1(end),'b*')
    plot(myi2,mxi2,'g-')
    plot(myi2(end),mxi2(end),'g*')

end
