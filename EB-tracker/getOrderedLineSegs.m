function [lineSegs,colMask]=getOrderedLineSegs

% make an example image
imL=512;
imW=512;

% here's one circle
centerYX=[100,200];
radius=35;

bwMask1=zeros(imL,imW);
bwMask1(centerYX(1),centerYX(2))=1;
bwMask1(bwdist(bwMask1)<radius)=1;
bwMask1=bwmorph(bwMask1,'remove');

% here's another circle to intersect with the first
centerYX=[130,240];
radius=50;

bwMask2=zeros(imL,imW);
bwMask2(centerYX(1),centerYX(2))=1;
bwMask2(bwdist(bwMask2)<radius)=1;
bwMask2=bwmorph(bwMask2,'remove');

% combine them
bwMask=bwMask1 | bwMask2;

% add a line
[myline,mycoords,outmat,X,Y]=bresenhamMatEx(bwMask1,[50 150; 300 350],0);
bwMask(sub2ind([imL imW],Y,X))=1;


% put a few breaks in the circles
[rows,cols]=find(bwMask);
pixInd=find(bwMask);

bwMask(rows(5),cols(5))=0;
rows(5)=[]; cols(5)=[];
pixInd(5)=[];

bwMask(rows(25),cols(25))=0;
rows(25)=[]; cols(25)=[];
pixInd(25)=[];

bwMask(rows(50),cols(50))=0;
rows(50)=[]; cols(50)=[];
pixInd(50)=[];

bwMask(rows(55),cols(55))=0;
rows(55)=[]; cols(55)=[];
pixInd(55)=[];

bwMask(rows(150),cols(150))=0;
rows(150)=[]; cols(150)=[];
pixInd(150)=[];

% end of example
imshow(bwMask)

bwMask=double(bwMask);
tempMask=bwMask; % use for plotting later

% get pix indices of BW image
[rows,cols]=find(bwMask==1);
pixInd=find(bwMask);

% neighbors are 1 or sqrt(2) pixels from each other if they are 4- or
% 8-connected, respectively
D=createDistanceMatrix([rows cols],[rows cols]);
D(D>sqrt(2))=0;

tooManyNeighborsIdx=find(sum(D>0,2)>2); % if 3 or more

bwMask(pixInd(tooManyNeighborsIdx))=.5; imshow(bwMask);


% take out those with too many neighbors, find D again
bwMask(pixInd(tooManyNeighborsIdx))=0;
[rows,cols]=find(bwMask);
pixInd=find(bwMask);
D=createDistanceMatrix([rows cols],[rows cols]);
D(D>sqrt(2))=0;


% possible endpoints for lines
segmentEnds=find(sum(D,2)<=sqrt(2));

segCount=1;

while ~isempty(segmentEnds)
    seg=segmentEnds(1);

    nextIdx=1;
    while ~isempty(nextIdx)

        cands=find(D(seg(end),:));
        nextIdx=setdiff(cands,seg);
        if ~isempty(nextIdx)
            seg=[seg nextIdx(1)];
        else
            segmentEnds=setdiff(segmentEnds,seg);
        end

    end

    lineSegs.pixList{segCount,1}=pixInd(seg);
    [y,x]=ind2sub([imL,imW],pixInd(seg));


    lineSegs.x{segCount,1}=x;
    lineSegs.y{segCount,1}=y;
    
    u=x(3:end)-x(1:end-2);
    v=y(3:end)-y(1:end-2);
    lineSegs.u{segCount,1}=zeros(size(x));
    lineSegs.v{segCount,1}=zeros(size(x));
    lineSegs.u{segCount,1}(2:end-1)=u;
    lineSegs.v{segCount,1}(2:end-1)=v;
    
    lineSegs.ori{segCount,1}=zeros(size(x));
    
    % % this gives 0=pi from top orientation
    %    lineSegs.ori{segCount,1}(2:end-1)=(pi()/2)-atan(-v./u);
    % % these lines give pi value when pointing up towards 0
    %     pointUpPix=find(u==0 & v~=0);
    %     lineSegs.ori{segCount,1}(pointUpPix)=pi();
    segCount=segCount+1;
end

nSegs=length(lineSegs.pixList);
tooShortIdx=[];
for i=1:nSegs
    tooShort=length(lineSegs.pixList{i,1})<3;
    if tooShort
    tooShortIdx=[tooShortIdx i];
    end
end
temp=lineSegs;
clear lineSegs

lineSegs.pixList=temp.pixList(setdiff(1:nSegs,tooShortIdx),1);
lineSegs.x=temp.x(setdiff(1:nSegs,tooShortIdx),1);
lineSegs.y=temp.y(setdiff(1:nSegs,tooShortIdx),1);
lineSegs.u=temp.u(setdiff(1:nSegs,tooShortIdx),1);
lineSegs.v=temp.v(setdiff(1:nSegs,tooShortIdx),1);
lineSegs.ori=temp.ori(setdiff(1:nSegs,tooShortIdx),1);


keptPix=cat(1,lineSegs.pixList{:});
tempMask(keptPix)=2;
tempMask=tempMask/max(tempMask(:));
figure(1); imshow(tempMask);

colMask=zeros(imL,imW,3);
r=zeros(imL,imW);
g=zeros(imL,imW);
b=zeros(imL,imW);

cmap=jet(length(lineSegs.pixList));
for i=1:length(lineSegs.pixList)
    r(lineSegs.pixList{i,1})=cmap(i,1);
    g(lineSegs.pixList{i,1})=cmap(i,2);
    b(lineSegs.pixList{i,1})=cmap(i,3);



    % cmap=jet(length(lineSegs(i,1).pixList));
    %
    % r(lineSegs(i,1).pixList)=cmap(:,1);
    % g(lineSegs(i,1).pixList)=cmap(:,2);
    % b(lineSegs(i,1).pixList)=cmap(:,3);

    % c=[rand rand rand];
    %
    % r(lineSegs(i,1).pixList)=c(1);
    % g(lineSegs(i,1).pixList)=c(2);
    % b(lineSegs(i,1).pixList)=c(2);

    colMask(:,:,1)=r;
    colMask(:,:,2)=g;
    colMask(:,:,3)=b;

    
end
figure(2);

imshow(colMask)
hold on
quiver(cat(1,lineSegs.x{:}),cat(1,lineSegs.y{:}),cat(1,lineSegs.u{:}),cat(1,lineSegs.v{:}),.25,'w');

oriMap=zeros(imL,imW,3);
r=zeros(imL,imW);
g=zeros(imL,imW);
b=zeros(imL,imW);

cmap2=jet(256);

oriValuesNorm=round(255.*cat(1,lineSegs.ori{:})./pi()+1);
oriValuesNorm(oriValuesNorm<1)=1; oriValuesNorm(oriValuesNorm>256)=256;

r(cat(1,lineSegs.pixList{:}))=cmap2(oriValuesNorm,1);
g(cat(1,lineSegs.pixList{:}))=cmap2(oriValuesNorm,2);
b(cat(1,lineSegs.pixList{:}))=cmap2(oriValuesNorm,3);

oriMap(:,:,1)=r;
oriMap(:,:,2)=g;
oriMap(:,:,3)=b;
figure(3); imshow(oriMap,[]);
disp('blah');





