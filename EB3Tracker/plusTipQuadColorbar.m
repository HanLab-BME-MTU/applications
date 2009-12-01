function [nPrctRYGB]=plusTipQuadColorbar(nPopRYGB)


nRows=size(nPopRYGB,1);
nPrctRYGB=zeros(nRows,4);
cBarAll=zeros(20*nRows,100,3);

cBCount=1;
for iRow=1:nRows
    popRYGB=nPopRYGB(iRow,:);

    if sum(popRYGB)==0
        cBCount=cBCount+20;
        continue
    end

    % pr is the rounded percent of red, yellow, green, and blue tracks there
    % are in the population
    percentsRYGB=round(100*(popRYGB./repmat(sum(popRYGB,2),[1 4])));

    % if percents don't add up to 100 exactly, subtract the difference from the
    % highest percentage
    d=sum(percentsRYGB,2);
    if d>100
        maxIdx=find(percentsRYGB==max(percentsRYGB,[],2),1);
        percentsRYGB(maxIdx)=percentsRYGB(maxIdx)-(d-100);
    elseif d<100
        minIdx=find(percentsRYGB==min(percentsRYGB,[],2),1);
        percentsRYGB(minIdx)=percentsRYGB(minIdx)+(100-d);
    end

    if sum(percentsRYGB,2)~=100
        error('Percents do not sum to 100')
    end

    cbarImg=zeros(20,100,3);

    % initialize black border
    b=ones(20,100);
    b(:,1)=0; b(:,end)=0; b(1,:)=0; b(end,:)=0;
    b=repmat(b,[1,1,3]);

    c=1;
    for i=1:4
        n=percentsRYGB(i);
        if n~=0
            switch i
                case 1
                    cbarImg(:,c:c+n-1,1)=1;
                case 2
                    cbarImg(:,c:c+n-1,1:2)=1;
                case 3
                    cbarImg(:,c:c+n-1,2)=1;
                case 4
                    cbarImg(:,c:c+n-1,3)=1;
            end
            c=c+n;
        end

    end



    nPrctRYGB(iRow,:)=percentsRYGB;
    cBarAll(cBCount:cBCount+20-1,:,:)=cbarImg.*b;

    cBCount=cBCount+20;

end

figure
imagesc(cBarAll);
axis equal
