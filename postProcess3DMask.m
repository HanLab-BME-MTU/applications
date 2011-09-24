function currMask = postProcess3DMask(currMask,p)


%if nargin < 2 || isempty(p)
    
p.MinVolume = 10;
p.ClosureRadius = 3;
p.FillHoles = 2;
p.NumObjects = 1;
p.FillDilateDiam = 5;
p.FuzzyFillThresh = 60;
p.CylinderFillDiam = 12;
p.CylinderFillHt = 30;
p.SuppressBorder = 6;

if p.SuppressBorder
    nSub = p.SuppressBorder;
    currMask([1:nSub end-nSub:end],:,:) = false;
    currMask(:,[1:nSub end-nSub:end],:) = false;
end


%end
closeBall = binarySphere(p.ClosureRadius);

currMask = bwareaopen(currMask,p.MinVolume);

currMask = imclose(currMask,closeBall);

labelMask = bwlabeln(currMask);



if p.NumObjects > 0
    rProp = regionprops(labelMask,'Area');
    [goodObj,iGood] = sort([rProp(:).Area],'descend'); %#ok<ASGLU>

    currMask = false(size(currMask));
    for i = 1:min(p.NumObjects,numel(iGood))
        currMask = currMask | (labelMask == iGood(i));
    end
end

if p.FillHoles

    if p.FillHoles > 1

        diSE = strel('disk',p.FillDilateDiam);

        %'Lenient' semi-2D hole filling                
        mFill1 = false(size(currMask));
        mFill2 = false(size(currMask));
        mFill3 = false(size(currMask));                
        for j = 1:size(currMask,1)
            %mFill1(j,:,:) = imfill(squeeze(currMask(j,:,:)),'holes');                    
            tmpDilated = imdilate(squeeze(currMask(j,:,:)),diSE);
            tmpHoles = imfill(tmpDilated,'holes') ~= tmpDilated;
            mFill1(j,:,:) = squeeze(currMask(j,:,:)) | imdilate(tmpHoles,diSE);        
        end
        for j = 1:size(currMask,2)
            %mFill2(:,j,:) = imfill(squeeze(currMask(:,j,:)),'holes');               
            tmpDilated = imdilate(squeeze(currMask(:,j,:)),diSE);
            tmpHoles = imfill(tmpDilated,'holes') ~= tmpDilated;
            mFill2(:,j,:) = squeeze(currMask(:,j,:)) | imdilate(tmpHoles,diSE);        
        end                    
        for j = 1:size(currMask,3)
            %mFill3(:,:,j) = imfill(currMask(:,:,j),'holes');                    
            tmpDilated = imdilate(squeeze(currMask(:,:,j)),diSE);
            tmpHoles = imfill(tmpDilated,'holes') ~= tmpDilated;
            mFill3(:,:,j) = squeeze(currMask(:,:,j)) | imdilate(tmpHoles,diSE);        
        end   
        %Only fill pixels which were filled in the specified number
        %of dimensions
        currMask = double(mFill1)+double(mFill2)+double(mFill3) >= p.FillHoles;                
        
        
        %TEMP - Cheating - the holes are almost always in the z direction, so we
        %bias the hole filling in this direction
        %currMask = double(mFill1)+double(mFill2)+2*double(mFill3) >= p.FillHoles;                

        %First do fuzzy hole filling
        maxDist = bwMaxDirectDist3D(currMask);
        maxDistSum = sum(maxDist,4);
        currMask = currMask | maxDistSum < p.FuzzyFillThresh;   
        
        %And cylinder filling
        cylDiamThresh = sum(maxDist<p.CylinderFillDiam,4);
        cylHtThresh = sum(maxDist<p.CylinderFillHt,4);
        
        %Fill in "cylindrical" areas smaller than the specified dimensions
        cylPoints = cylDiamThresh  >= 4 & cylHtThresh >= 5;
        
        currMask = currMask | cylPoints;
        
        %Follow this with another round of closure
        currMask = imclose(currMask,closeBall);                        

    end
    
    %Do strict, full-3D hole-filling which includes the image border as a
    %boundary
    CC = bwconncomp(~currMask,6);
    %Keep only the largest background object
    if CC.NumObjects > 1
        [~,iBySize] = sort(cellfun(@numel,CC.PixelIdxList));
        for j = 1:CC.NumObjects-1
            currMask(CC.PixelIdxList{iBySize(j)}) = true;
        end
    end
    
    
end

