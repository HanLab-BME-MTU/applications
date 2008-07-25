function eb1SpotDetector(runInfo,nFrames,debug)

% runInfo: structure containing image and analysis directories
% nFrames: number of frames on which to detect; []=all (default)
% debug: 0 to detect spots on all images; 1 to test one image

warningState=warning;
warning('off','MATLAB:divideByZero')

if nargin<3
    error('eb1SpotDetector: need 3 input arguments')
end

% check runInfo format and directory names between file systems
if ~isfield(runInfo,'imDir') || ~isfield(runInfo,'anDir')
    error('eb1SpotDetector: runInfo should contain fields imDir and anDir');
else
    [runInfo.anDir]=formatPath(runInfo.anDir);
    [runInfo.imDir]=formatPath(runInfo.imDir);
end

% get roi edge pixels and make region outside mask NaN
if ~isfield(runInfo,'roiMask')
    roiMask=1;
    edgePix=[];
else
    polyEdge=bwmorph(runInfo.roiMask,'remove');
    edgePix=find(polyEdge);
    roiMask=swapMaskValues(double(runInfo.roiMask),0,NaN);
end

% make spot directory if it doesn't exist from batch
spotDir=[runInfo.anDir filesep 'spot'];
if ~isdir(spotDir)
    mkdir(spotDir);
end

% count number of images in image directory
[listOfImages] = searchFiles('.tif',[],runInfo.imDir,0);
nImTot=size(listOfImages,1);

if debug==1 % let user select image to test
    currentDir=pwd;
    cd(runInfo.imDir)
    [fileName,dirName] = uigetfile('*.tif','Choose a .tif file');
    cd(currentDir)

    for i=1:nImTot
        if strcmp(listOfImages{i,1},fileName)
            break
        end
    end
    m=i;
    n=i;

elseif debug==0 % do detection on all images
    m=1;
    if ~isempty(nFrames) && nFrames>0 && nFrames<nImTot
        n=nFrames;
    else
        n=nImTot;
    end
end


% sigma for gauss filtering
gaussSigma=4;


% string for sigma in gauss filtered
s1=length(num2str(max(gaussSigma)));
strg1=sprintf('%%.%dd',s1);

% string for number of files
s2=length(num2str(length(m:n)));
strg2=sprintf('%%.%dd',s2);

% string for percent max cutoff
s3=length(num2str(2));
strg3=sprintf('%%.%dd',s3);

% get image dimensions
[imL,imW]=size(double(imread([char(listOfImages(1,2)) filesep char(listOfImages(1,1))])));

% agg=zeros(imL,imW,length(gaussSigma));
% for sig=1:length(gaussSigma)
%         fileNameIm=[char(listOfImages(m,2)) filesep char(listOfImages(m,1))];
%         img=double(imread(fileNameIm));
%
%         lowPass=Gauss2DBorder(img,1);
%         highPass=Gauss2DBorder(img,gaussSigma(sig));
%         filterDiff=lowPass-highPass;
%
%         agg(:,:,sig)=filterDiff;
% end
% figure(1); imshow(sum(agg,3),[])
% figure(2); imshow(mean(agg,3),[])


% initialize structure to store info for tracking
movieInfo(n-m+1,1).xCoord=0;
movieInfo(n-m+1,1).yCoord=0;
movieInfo(n-m+1,1).amp=0;
movieInfo(n-m+1,1).int=0;

for i=m:n % loop thru frames

    if debug==1
        agg=zeros(imL,imW,length(gaussSigma));
        minMax=zeros(length(gaussSigma),2);
    end

    for sig=1:length(gaussSigma) % loop thru sigma sizes

        indxStr1=sprintf(strg1,gaussSigma(sig));

        fileNameIm=[char(listOfImages(i,2)) filesep char(listOfImages(i,1))];
        img=double(imread(fileNameIm));

        lowPass=Gauss2DBorder(img,1);
        highPass=Gauss2DBorder(img,gaussSigma(sig));
        filterDiff=roiMask.*(lowPass-highPass);
        %filterDiff=roiMask.*(lowPass-meanImg);
        
        %if debug==1
        %         minMax(sig,:)=[min(filterDiff(:)) max(filterDiff(:))];
        %         maxAbs=max(abs(minMax(sig,:)));
        %         normFilterDiff=filterDiff./maxAbs;
        %
        %         RGB=zeros(imL,imW,3);
        %         r=zeros(imL,imW); r(normFilterDiff>=0)=normFilterDiff(normFilterDiff>=0);
        %         g=zeros(imL,imW); g(normFilterDiff<=0)=abs(normFilterDiff(normFilterDiff<=0));
        %         RGB(:,:,1)=r;
        %         RGB(:,:,2)=g;
        %         RGB=round(255.*(RGB.^gamma));

        %         imwrite(RGB,[spotDir filesep 'filterDiff_RG_sigma_' num2str(sig) '.tif']);
        %         agg(:,:,sig)=filterDiff;

        %     agg=reshape(agg,[imL*imW,length(gaussSigma)]);
        %     edges=round(linspace(min(minMax(:)),max(minMax(:)),100)); % 100
        %     [n1,xout1]=histc(agg(:,1),edges);
        %     [n2,xout2]=histc(agg(:,2),edges);
        %     [n3,xout3]=histc(agg(:,3),edges);
        %     [n4,xout4]=histc(agg(:,4),edges);
        %     [n5,xout5]=histc(agg(:,5),edges);
        %     [n6,xout6]=histc(agg(:,6),edges);
        %     Y=[n1'; n2'; n3'; n4'; n5'; n6';];

        %     % histogram
        %     bar(edges,Y','group')
        %     legend(num2str(gaussSigma'))
        %     title(['frame ' num2str(i) ': Gauss2D(img,1) - Gauss(img,sigma) intensity histogram'])

        %     % cumulative histogram
        %     figure(2); cY=cumsum(Y');
        %     bar(edges,cY,'group') %whatever, fix this tomorrow...
        %     legend(num2str(gaussSigma'))
        %     title(['frame ' num2str(i) ': Gauss2D(img,1) - Gauss(img,sigma) intensity cumulative histogram'])
        %
        %end

        maxFiltDiff=nanmax(filterDiff(:));
        minFiltDiff=nanmin(filterDiff(:));

        % assume lowest value we care about is x% of max value of
        % filterDiff, where x is the abs(most negative value)
        %lowCutPercent=round(100*(abs(minFiltDiff)/maxFiltDiff));
         lowCutPercent=33;
        % compare features in z-slices startest from the highest one
        count=1;
        for p=99:-1:lowCutPercent % percent of max intensity of filterDiff

            if count==1
                slice1=filterDiff>.01*p*maxFiltDiff; % top slice
                [featMap1,numFeats1]=bwlabel(slice1);
                featProp1=regionprops(featMap1,'PixelIdxList','Area');
            else
                featProp1=featProp2;
                numFeats1=numFeats2;
            end
            slice2=filterDiff>.01*(p-1)*maxFiltDiff; % next slice down
            imshow(slice2,[])
            [featMap2,numFeats2]=bwlabel(slice2);
            featProp2=regionprops(featMap2,'PixelIdxList','Area','Extrema');

            % groupIntersect looks for pixls shared by features in slice1 and
            % slice2. since pixels can only belong to one feature, as we move
            % down the z-slices, we will find that features in the lower slice
            % have intersections with more than one feature from the upper
            % slice. we want to kick out the merged feature and replace it with
            % the smaller individual features from the upper slice.
            groupIntersect=zeros(numFeats1,numFeats2);
            for rows=1:numFeats1
                for cols=1:numFeats2
                    groupIntersect(rows,cols)=~isempty(intersect(featProp1(rows,1).PixelIdxList,featProp2(cols,1).PixelIdxList));
                end
            end

            % index of features in slice1 that have more than one feature from
            % slice2 contributing
            mergeIdx=find(sum(groupIntersect,1)>1);

            for idx=1:length(mergeIdx)
                % replace pixel area with NaN for a merged feature in slice2
                featProp2(mergeIdx(idx)).Area=nan;
                % get indices for the features from slice1 that contributed to
                % the one we just kicked out and add these to the end of the
                % slice2 feature list
                rowIdx=find(groupIntersect(:,mergeIdx(idx)));
                for rows=1:length(rowIdx)
                    numFeats2=numFeats2+1; % increase feature number since added
                    featProp2(numFeats2)=featProp1(rowIdx(rows));
                end
            end

            count=count+1;
        end

        %smallFeats=find([featProp2(:,1).Area]<=4);

        % the good features have area (weren't a merged feature; ~NaN) and we
        % assume the area > 4 pixels
        goodFeats=find(~isnan([featProp2(:,1).Area]))';% & [featProp2(:,1).Area]>4)';
        featureList=deal(featProp2(goodFeats));
        
        featureMap=zeros(imL,imW);
        numFeats=length(goodFeats);
        featureIntensity=zeros(numFeats,1);
        for iFeat=1:length(goodFeats)
            % make new labeled map 
            featureMap(featureList(iFeat,1).PixelIdxList)=iFeat;
            
            % integrate intensity for each feature
            featureIntensity(iFeat)=sum(img(featureList(iFeat,1).PixelIdxList));
        end
        %get feature properties
        featPropFinal=regionprops(featureMap,'PixelIdxList','Area','Centroid','Extrema'); %,'Orientation');

        % centroid coordinates
        yCoord=zeros(numFeats,2); xCoord=zeros(numFeats,2);
        pos=reshape([featPropFinal.Centroid],2,[])'; 
        yCoord(:,1)=pos(:,2); xCoord(:,1)=pos(:,1);
        
        % area
        featArea=[featPropFinal.Area]';
        
        amp=ones(numFeats,2); % assume area std is 1??
        amp(:,1)=featArea;
        
        
        % orientation in radians
%         ori=[featPropFinal.Orientation]'.*pi/180;
%         u=5*-cos(ori);
%         v=5*sin(ori);

        % use extrema to draw polygon around spots - here we get
        % coordinates for polygon
        outline=[featPropFinal.Extrema]; outline=[outline; outline(1,:)];
        outlineX=outline(:,1:2:size(outline,2));
        outlineY=outline(:,2:2:size(outline,2));

        
        indxStr2=sprintf(strg2,i); %frame
        indxStr3=sprintf(strg3,lowCutPercent); %percentMaximum
 
        % plot spot outlines and centroid on image
        figure(1); 
        im2show=lowPass./max(abs(lowPass(:))); % image for overlay
        im2show(edgePix)=1; % ROI in white
        imshow(im2show,[]);
        hold on
        scatter(xCoord(:,1),yCoord(:,1),'c.'); % plot centroid in cyan
        %quiver(xCoord(:,1),yCoord(:,1),u,v,0,'r'); % plot directionality in red
        plot(outlineX,outlineY); % plot spot outlines in random colors
 
        saveas(gcf,[spotDir filesep 'overlay_f' indxStr2 '_p' indxStr3])
        saveas(gcf,[spotDir filesep 'overlay_f' indxStr2 '_p' indxStr3 '.tif']);
        
        
        % make color image of spots
        RGB=label2rgb(featureMap, 'jet', 'k','shuffle');
        
        % save image of colored spots
        figure(2); 
        imshow(RGB); 
        hold on;
        %scatter(xCoord(:,1),yCoord(:,1),'c.');
        %quiver(xCoord(:,1),yCoord(:,1),u,v,0,'r');
        %plot(outlineX,outlineY);
        saveas(gcf,[spotDir filesep 'spots_f' indxStr2 '_p' indxStr3 '.tif']);
        
        
        % make structure compatible with Khuloud's tracker
        movieInfo(i,1).xCoord=xCoord;
        movieInfo(i,1).yCoord=yCoord;
        movieInfo(i,1).amp=amp;
        movieInfo(i,1).int=featureIntensity;
        save([spotDir filesep 'movieInfo'],'movieInfo');
        
        % save a structure per frame
        %save([spotDir filesep 'featureList_f' indxStr2],'featPropFinal');
    end

end
warning(warningState);



% go thru these features and look at area before/after
%         for idx=1:length(mergeIdx)
%             feat2merge=find(groupIntersect(:,mergeIdx(idx)));
%             feats1AreaSum=0;
%             for f=1:length(feat2merge)
%                 feats1AreaSum=feats1AreaSum+featProp1(feat2merge(f),1).Area;
%                 feats2Area=featProp2(mergeIdx(idx),1).Area;
%             end
%         end











% [n,xout]=hist(filterDiff(:),100);
%
% for j=1:100
%     imshow(filterDiff>xout(97))
%     bwSpotMap=filterDiff>xout(end-j+1);
%     [featMap,numFeats]=bwlabel(bwSpotMap,4);
%
%     RGB=label2rgb(featMap, 'jet', 'k','shuffle');
%     h=gcf;
%     figure(h); imshow(RGB);
%     imwrite(RGB,[spotDir filesep 'junk' num2str(j) '.tif']);
% end
%
% [cutoffInd,cutoffValue]=cutFirstHistMode(filterDiff,1);
% drawnow
% % coef = 4 Katsu; coef = 1 Claudio; coef = 1 Lisa_xju103_r11;
% bwSpotMap=filterDiff>cutoffValue*coef; % REMOVE THE NOISE FEATURES %no 3
%
% [featMap,numFeats]=bwlabel(bwSpotMap);
% RGB=label2rgb(featMap, 'jet', 'k','shuffle');
% h=gcf;
% figure(h+1); imshow(RGB);
%
% warningState=warning;
% %    intwarning off
% feats=regionprops(featMap,'Centroid','Eccentricity','Orientation','MajorAxisLength');
% warning(warningState);
%
% %     % Initialize 'feats' structure
% %     feats(numFeats)=struct(...
% %         'pos',[0 0],...
% %         'ecc',0,...
% %         'ori',0,...
% %         'len',0);
%
% figure(4); imshow(img,[]); hold on
% for j=1:length(feats)
%     e1 = [-cos(feats(j).Orientation*pi/180) sin(feats(j).Orientation*pi/180) 0];
%     e2 = [sin(feats(j).Orientation*pi/180) cos(feats(j).Orientation*pi/180) 0];
%     e3 = [0 0 1];
%     Ori = [feats(j).Centroid  0];
%     v1 = [-10 10];
%     v2 = [-5 5];
%     v3 = [0 0];
%     [xGrid,yGrid]=arbitraryGrid(e1,e2,e3,Ori,v1,v2,v3);
%
%     imshow(featMap==j); hold on
%     scatter(xGrid(:),yGrid(:),'r+');
%
%     Crop(:,:,j) = interp2(img,xGrid,yGrid);
%     %         Crop(:,:,j) = interp2(img,xGrid,yGrid,'*linear');
%
%     e1 = [];e2 = [];e3 = []; Ori = []; v1 = []; v2 = []; xGrid = []; yGrid = [];
% end
%
% Cm = nanmean(Crop,3); % MEAN/REPRESENTATIVE EB1 CROP
% Crop(isnan(Crop))=0;% border effect - some NaN
% Cm1 = bwlabel(Cm);
% statsC = regionprops(Cm1,'all');
%
% sC = size(Crop);
% Cm3d = repmat(Cm,[1,1,size(Crop,3)]);
% dC = Crop - Cm3d;
% sqC = dC.^2;
% ssqC = squeeze(sum(sum(sqC,1),2)); %LIST OF DIFFERENCES AFTER SUBTRACTION
%
% B = Cm(:);
% A = ones(length(B),2);
%
% for m = 1:size(Crop,3)
%     CR = Crop(:,:,m);
%     A(:,2) = CR(:);
%     goodRows = find(A(:,2) ~= 0 & isfinite(B));
%     XX = lscov(A(goodRows,:),B(goodRows));
%     RES = B(goodRows) - A(goodRows,:)*XX;
%     OUT(m,:) = [mean(RES(:).^2),XX'];
% end
%
% [Ind,V]=cutFirstHistMode(OUT(:,1),0);% switch to 1 to see HIST
%
% goodFeats = find(OUT(:,1)<V); % SPOTS WHICH FIT WELL WITH THE MEAN EB1 SPOT
%
% featNames = fieldnames(feats);
% for field = 1:length(featNames)
%     feats.(featNames{field}) = feats.(featNames{field})(goodFeats,:);
% end
%
% if debug == 1
%     alowPass = 5;
%     If=Gauss2D(img,1);
%     figure, imshow(If(1+alowPass:end-alowPass,1+alowPass:end-alowPass),[]);%bwSpotMap
%     hold on
%     for img = 1:length(feats.ori)
%         h = quiver(feats.pos(img,1)-alowPass,feats.pos(img,2)-alowPass,-cos(feats.ori(img)*pi/180),sin(feats.ori(img)*pi/180),3,'r');
%         set(h,'LineWidth',2)
%     end
% elseif debug == 0
%     save([dirName(1:end-7),'cands',filesep,'feats',indxStr],'feats')
%     clear goodFeats
%     clear OUT
%     clear V
%     clear Crop
% end



