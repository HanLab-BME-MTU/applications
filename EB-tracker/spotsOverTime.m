function spotsOverTime(runInfo,nFrames2Plot)

% get list of featureList files
temp=load([runInfo.anDir filesep 'spot' filesep 'movieInfo']);
movieInfo=temp.movieInfo;
nDetectedFrames=size(movieInfo,1);


% string for number of files
s=length(num2str(nDetectedFrames));
strg=sprintf('%%.%dd',s);


% get list of images in integrated image directory
intImDir=[runInfo.anDir filesep 'spot' filesep 'intIm'];
[listOfImages] = searchFiles('meanImg',[],intImDir,0);
nImTot=size(listOfImages,1);


cMap=jet(nFrames2Plot);
for i=1:nDetectedFrames-nFrames2Plot+1
    
   img=double(imread([char(listOfImages(i,2)) filesep char(listOfImages(i,1))]));
   %figure(1); imshow(img,[]);
   hold on
   
   for j=i:i+nFrames2Plot-1
       % load feature list
       scatter(movieInfo(j,1).xCoord(:,1),movieInfo(j,1).yCoord(:,1),[],cMap(j-i+1,:)); % plot centroid 
       
   end
   
   indxStr=sprintf(strg,i);
   saveas(gcf,[intImDir filesep 'centroid_overlay_f' indxStr]);
   saveas(gcf,[intImDir filesep 'centroid_overlay_f' indxStr '.tif']);
   
end





