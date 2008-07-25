function getCircle

imL=512;
imW=512;

centerYX=[100,200];
radius=35;


close

cMask=zeros(imL,imW);
cMask(centerYX(1),centerYX(2))=1;
cMask(bwdist(cMask)<radius)=1;
cMask=bwmorph(cMask,'remove');


[rows,cols]=find(cMask);
pixInd=find(cMask);

cMask(rows(5),cols(5))=0;
rows(5)=[]; cols(5)=[];
pixInd(5)=[];

cMask(rows(50),cols(50))=0;
rows(50)=[]; cols(50)=[];
pixInd(50)=[];

cMask(rows(25),cols(25))=0;
rows(25)=[]; cols(25)=[];
pixInd(25)=[];

cMask(rows(55),cols(55))=0;
rows(55)=[]; cols(55)=[];
pixInd(55)=[];


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
    
    lineSegs(segCount,1).pixList=pixInd(seg);
    
    segCount=segCount+1;
end


figure(1); imshow(cMask);

colMask=zeros(imL,imW,3);
r=zeros(imL,imW); 
g=zeros(imL,imW); 
b=zeros(imL,imW); 

for i=1:length(lineSegs)

% cmap=jet(length(lineSegs(i,1).pixList));
% 
% r(lineSegs(i,1).pixList)=cmap(:,1);
% g(lineSegs(i,1).pixList)=cmap(:,2);
% b(lineSegs(i,1).pixList)=cmap(:,3);
c=[rand rand rand];

r(lineSegs(i,1).pixList)=c(1);
g(lineSegs(i,1).pixList)=c(2);
b(lineSegs(i,1).pixList)=c(2);

colMask(:,:,1)=r;
colMask(:,:,2)=g;
colMask(:,:,3)=b;

end
figure(2); imshow(colMask)
