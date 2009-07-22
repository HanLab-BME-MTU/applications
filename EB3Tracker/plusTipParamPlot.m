function [thresh1,thresh2]=plusTipParamPlot(param1,cutoff1,param2,cutoff2)

if nargin<1
    param1='growthSpeed';
    cutoff1=50;
    param2='growthLifetime';
    cutoff2=50;
end

dirName=uigetdir(pwd,'Please select project directory');
temp=load([dirName filesep 'meta' filesep 'projData.mat']);
projData=temp.projData;

% format image/analysis directory paths
imDir=formatPath(projData.imDir);
anDir=formatPath(projData.anDir);

% get first image from imDir
[listOfImages] = searchFiles('.tif',[],imDir,0);
img = double(imread([listOfImages{1,2} filesep listOfImages{1,1}]));

% extract track info (aggregated)
allData=abs(projData.nTrack_sF_eF_vMicPerMin_trackType_lifetime_totalDispPix);

% get matrices containing coordinates for all subtracks
trackType=allData(:,5);
[xMat,yMat]=plusTipGetSubtrackCoords(projData,[]);


for i=1:2
    % assign data based on user input
    subtrackIdx=find(trackType==1);
    switch eval(['param' num2str(i)])
        case 'growthSpeed'
            temp=allData(subtrackIdx,4);
        case 'growthLifetime'
            temp=allData(subtrackIdx,6);
        otherwise
            temp='This parameter is not supported';
    end

    if i==1
        data1=temp;
    else
        data2=temp;
    end

end

% figure; hist(data1,25)
% title(param1);
%
% figure; hist(data2,25)
% title(param2);

colorMap=['b','g','y','r'];

figure(2);
imagesc(img); colormap gray;

thresh1=prctile(data1,cutoff1);
thresh2=prctile(data2,cutoff2);

for i=1:4
    switch i
        case 1 % fast and long
            idx=find(data1>thresh1 & data2>thresh2);
        case 2 % slow and long
            idx=find(data1<=thresh1 & data2>thresh2);
        case 3 % fast and short
            idx=find(data1>thresh1 & data2<=thresh2);
        case 4 % slow and short
            idx=find(data1<=thresh1 & data2<=thresh2);     
    end

    figure(1); hold on
    scatter(data1(idx),data2(idx),[],colorMap(i),'.');
    xlabel(param1)
    ylabel(param2)

    figure(2); hold on
    plot(xMat(subtrackIdx(idx),:)',yMat(subtrackIdx(idx),:)',colorMap(i))

end



