function projImages=project1D(MD,dynPoligon,varargin)
% tracks are in original FoF and projected  in the manifold
% described by polygon.
% Poligon described the projected space (without dilated boundaries), It is composed of Tracks.
% in 1D, it is only composed of two tracks.

ip = inputParser;
ip.CaseSensitive = false;
ip.KeepUnmatched = true;
ip.addOptional('tracks',[]);
ip.parse(varargin{:});
p=ip.Results;

tracks=p.tracks;
%% Test on a single kinetochore bundle

% outputDirBundle=[MD.outputDirectory_ filesep 'Kin' filesep 'bundles'];
% tmp=load([outputDirBundle filesep 'kin-MT-bundle.mat'],'kinTracks');
% kinTracksBundle=tmp.kinTracks;
%%

showDebugGraphics=0;
cubeHalfWidth=20;


%% Cube that contain the moving polygon
minX=MD.getDimensions('X')+1;
minY=MD.getDimensions('Y')+1;
minZ=MD.getDimensions('Z')+1;
maxX=0;
maxY=0;
maxZ=0;
for iP=1:length(dynPoligon)
    minX=floor(min(min(dynPoligon(iP).x),minX));
    minY=floor(min(min(dynPoligon(iP).y),minY));
    minZ=floor(min(min(dynPoligon(iP).z),minZ));
    maxX=ceil(max(max(dynPoligon(iP).x),maxX));
    maxY=ceil(max(max(dynPoligon(iP).y),maxY));
    maxZ=ceil(max(max(dynPoligon(iP).z),maxZ));    
end
%%

outputDirSlices1=[MD.outputDirectory_ filesep '1DProjection' ];
outputDirDemo=[MD.outputDirectory_ filesep '1DProjection' filesep 'volDemo' ];

system(['mkdir  ' outputDirSlices1]);
system(['mkdir  ' outputDirDemo]);
for fIdx=1:MD.nFrames_
    vol=MD.getChannel(1).loadStack(fIdx);
    kinvol=MD.getChannel(2).loadStack(fIdx);
    % Collect relative frameIdx
    pIndices=nan(1,length(dynPoligon));
    for polIdx=1:length(dynPoligon)
        F=dynPoligon(polIdx).f;
        pIdx=find(F==fIdx);
        if isempty(pIdx)
            if(fIdx>max(F))   pIdx=max(F);  else   pIdx=min(F); end;
        end
        pIndices(polIdx)=pIdx;
    end;
    PCurrent=[dynPoligon(1).x(pIndices(1)) dynPoligon(1).y(pIndices(1)) dynPoligon(1).z(pIndices(1))];
    KCurrent=[dynPoligon(2).x(pIndices(2)) dynPoligon(2).y(pIndices(2)) dynPoligon(2).z(pIndices(2))];

    % Building mask for both channel on the whole volume
    mask=zeros(size(vol));
    sampling=100;
    xSeg=round(linspace(PCurrent(1),KCurrent(1),sampling));
    ySeg=round(linspace(PCurrent(2),KCurrent(2),sampling));
    zSeg=round(linspace(PCurrent(3)*MD.pixelSize_/MD.pixelSizeZ_,KCurrent(3)*MD.pixelSize_/MD.pixelSizeZ_,sampling));
    indx=sub2ind(size(mask),ySeg,xSeg,zSeg);

    mask(indx)=1;
    mask=imdilate(mask,ones(cubeHalfWidth,cubeHalfWidth,round(cubeHalfWidth*MD.pixelSize_/MD.pixelSizeZ_)));

    maskedVol=vol;
    maskedVol(~mask)=0;
    outputDir=[outputDirDemo filesep  'mask'];mkdir(outputDir);
    imwrite(uint8(255*mat2gray(max(maskedVol,[],3))),[outputDir filesep  'mask_' num2str(fIdx,'%04d') '.png']) % Names the file stitched001..., in /stitched/

    maskedKin=kinvol;
    maskedKin(~mask)=0;
    outputDir=[outputDirDemo filesep  'maskKin_'];mkdir(outputDir);
    imwrite(uint8(255*mat2gray(max(maskedKin,[],3))),[outputDir filesep  'maskKin_' num2str(fIdx,'%04d') '.png']) % Names the file stitched001..., in /stitched/


    % If needed the map must rotated before cropped (efficiency)
    % Rotation depends on the FrameOfRef associated to the tracks the compose the dynanimc polygon
    % Cropping area according to the polygon OVER TIME plus added vizualiation margin 
    % Rotation will use imwarp
    % Can we use imwar for cropping too ? 
    
    %% testing imwarp to crop the image
    % cropping zine
    maxXBorder=min(maxX+cubeHalfWidth,MD.getDimensions('X'));
    maxYBorder=min(maxY+cubeHalfWidth,MD.getDimensions('Y'));
    maxZBorder=min(maxZ+cubeHalfWidth,MD.getDimensions('Z'));
    minXBorder=max(minX-cubeHalfWidth,1);
    minYBorder=max(minY-cubeHalfWidth,1);
    minZBorder=max(minZ-cubeHalfWidth,1);    

    tform=affine3d();
    
    minZ=minZ*MD.pixelSize_/MD.pixelSizeZ_;
    maxZ=maxZ*MD.pixelSize_/MD.pixelSizeZ_;

    outputRef=imref3d([maxXBorder-minXBorder maxYBorder-minYBorder MD.getDimensions('Z')], [minXBorder maxXBorder],[minYBorder maxYBorder],[1 MD.getDimensions('Z')]);
    warpedVol=imwarp(vol,tform,'OutputView',outputRef);
    warpedMaskedVol=imwarp(maskedVol,tform,'OutputView',outputRef);
    warpedKinVol=imwarp(kinvol,tform,'OutputView',outputRef);
    warpedMaskedKinVol=imwarp(maskedKin,tform,'OutputView',outputRef);


%     figure
%     subplot(1,3,1);
%     imshow(mat2gray(max(vol,[],3)));
%     subplot(1,3,2);
%     imshow(mat2gray(max(warpedVol,[],3)));
%     subplot(1,3,3);
%     imshow(mat2gray(max(maskedVol,[],3)));
    imwrite(uint8(255*mat2gray(max(vol,[],3))),[outputDirDemo filesep  'vol_' num2str(fIdx,'%04d') '.png']) % Names the file stitched001..., in /stitched/
    imwrite(uint8(255*mat2gray(max(warpedVol,[],3))),[outputDirDemo filesep  'warpedVol_' num2str(fIdx,'%04d') '.png']) % Names the file stitched001..., in /stitched/
    imwrite(uint8(255*mat2gray(max(warpedMaskedVol,[],3))),[outputDirDemo filesep  'warpedMaskedVol_' num2str(fIdx,'%04d') '.png']) % Names the file stitched001..., in /stitched/
    imwrite(uint8(255*mat2gray(max(warpedKinVol,[],3))),[outputDirDemo filesep  'warpedKinVol_' num2str(fIdx,'%04d') '.png']) % Names the file stitched001..., in /stitched/
    imwrite(uint8(255*mat2gray(max(warpedMaskedKinVol,[],3))),[outputDirDemo filesep  'warpedKinMaskedVol_' num2str(fIdx,'%04d') '.png']) % Names the file stitched001..., in /stitched/

%     
    %%
    % Cropping mask to area of interest
    maskcrop=maskedVol;
    nullMaskXY=(squeeze(any(maskcrop,3)));
    YNull=~(squeeze(any(any(mask,3),2)));
    XNull=~(squeeze(any(any(mask,3),1)));
    ZNull=~(squeeze(any(any(mask,1),2)));
    maskcrop(:,:,ZNull)=[];
    maskcrop(YNull,:,:)=[];
    maskcrop(:,XNull,:)=[];
    maskedKin(:,:,ZNull)=[];
    maskedKin(YNull,:,:)=[];
    maskedKin(:,XNull,:)=[];

    if(showDebugGraphics)
        subplot(2,2,1);
        imshow(nullMaskXY)
        subplot(2,2,2);
        plot(YNull);
        subplot(2,2,3);
        plot(YNull);
        subplot(2,2,4);
        plot(ZNull);
        nullMaskXY=(squeeze(any(maskcrop,3)));
        subplot(2,2,1);
        imshow(nullMaskXY)
    end

    capturedEB3OrigRef=tracks;

    [maxXY,maxZY,maxZX,~]=computeMIPs(maskcrop,216/MD.pixelSize_,3*min(vol(:)),0.9*max(maskedVol(:)));
    [maxXYKin,maxZYKin,maxZXKin,~]=computeMIPs(maskedKin,216/MD.pixelSize_,0.6*max(maskedKin(:)),max(maskedKin(:)));

    
    
    
    %
    % Relative scaling

    maxMIPSize=400;
    [sX,sY,sZ]=size(maskcrop);
    sZ=sZ*MD.pixelSizeZ_/MD.pixelSize_;
    resizeScale=maxMIPSize/max([sX,sY,sZ]);


    rMaxXY=imresize(maxXY,resizeScale);
    rmaxXYKin=imresize(maxXYKin,resizeScale);
    mipSize=size(rMaxXY);
%     mipSize=[400 400];
%     rMaxXY=imresize(maxXY,mipSize);
%     rmaxXYKin=imresize(maxXYKin,mipSize);

    RGBThree=repmat(rMaxXY,1,1,3);
    % channel fusion
    RGBThree(:,:,1)=max(rMaxXY,rmaxXYKin);

    myColormap=[[0 0 255];
    [0 255 00]];


    tracksColors=uint8(myColormap);
    if(~isempty(tracks))
      tracksXY=trackBinaryOverlay(RGBThree,[find(~XNull,1) find(~XNull,1,'last')],[find(~YNull,1) find(~YNull,1,'last')],capturedEB3OrigRef,fIdx,ones(1,length(capturedEB3OrigRef)),tracksColors);
    else
      tracksXY=RGBThree;
    end

    if (size(tracksXY,2)<maxMIPSize)
        tracksXY=padarray(tracksXY,[0 maxMIPSize-size(tracksXY,2)],'post');
    end
    if (size(tracksXY,1)<maxMIPSize)
        tracksXY=padarray(tracksXY,[maxMIPSize-size(tracksXY,1) 0],'post');
    end


    rMax=imresize(maxZY,resizeScale);
    rmaxKin=imresize(maxZYKin,resizeScale);

    RGBThree=repmat(rMax,1,1,3);
    RGBThree(:,:,1)=max(rMax,rmaxKin);

    if(~isempty(tracks))
      capturedEB3ZY=capturedEB3OrigRef.copy();
      for ebIdx=1:length(capturedEB3ZY)
        capturedEB3ZY(ebIdx).x=capturedEB3OrigRef(ebIdx).z*MD.pixelSize_/MD.pixelSizeZ_;
      end
      tracksColors=uint8(myColormap);
      tracksZY=trackBinaryOverlay(RGBThree,[find(~ZNull,1) find(~ZNull,1,'last')],[find(~YNull,1) find(~YNull,1,'last')],capturedEB3ZY,fIdx,ones(1,length(capturedEB3OrigRef)),tracksColors);
    else
      tracksZY=RGBThree;
    end;    

    if (size(tracksZY,2)<maxMIPSize)
      tracksZY=padarray(tracksZY,[0 maxMIPSize-size(tracksZY,2)],'post');
    end;
    if (size(tracksZY,1)<maxMIPSize)
      tracksZY=padarray(tracksZY,[maxMIPSize-size(tracksZY,1)],'post');
    end;

    tracksZY=permute(tracksZY,[2 1 3]);
    stripeSize=4;
    threeTop = [tracksXY, zeros(size(tracksXY,1), stripeSize,3), zeros(size(tracksXY,1), size(tracksZY,2),3)];


    rMax=imresize(maxZX,resizeScale);
    rmaxKin=imresize(maxZXKin,resizeScale);

    RGBThree=repmat(rMax,1,1,3);
    RGBThree(:,:,1)=max(rMax,rmaxKin);

    if(~isempty(tracks))
      capturedEB3ZX=capturedEB3ZY.copy();
      for ebIdx=1:length(capturedEB3ZX)
        capturedEB3ZX(ebIdx).y=capturedEB3OrigRef(ebIdx).x;
      end
      tracksColors=uint8(myColormap);
      tracksZX=trackBinaryOverlay(RGBThree,[find(~ZNull,1) find(~ZNull,1,'last')],[find(~XNull,1) find(~XNull,1,'last')],capturedEB3ZX,fIdx,ones(1,length(capturedEB3OrigRef)),tracksColors);
    else
      tracksZX=RGBThree;
    end;

    if (size(tracksZX,2)<maxMIPSize)
      tracksZX=padarray(tracksZX,[0 maxMIPSize-size(tracksZX,2)],'post');
    end;
    if (size(tracksZX,1)<maxMIPSize)
      tracksZX=padarray(tracksZX,[maxMIPSize-size(tracksZX,1)],'post');
    end;


    tracksZX=permute(tracksZX,[2 1 3]);
    threeBottom = [tracksZX, 0*ones(size(tracksZX,1),+stripeSize,3),tracksZY];
    three = [threeTop; ones(stripeSize, size(tracksXY,2)+size(tracksZY,2)+stripeSize,3); threeBottom];
    %%
    imwrite(three,[outputDirSlices1 filesep 'frame_nb' num2str(fIdx,'%04d') '.png']);
end

video = VideoWriter([outputDirSlices1  '.avi']);
video.FrameRate = 5;  % Default 30
video.Quality = 100;    % Default 75

open(video)
for frameIdx=1:MD.nFrames_
    % save the maximum intensity projections
%     for poleIdx=poleIndices
%         outputDirSlices1=[MD.outputDirectory_ filesep 'maskSliceTracksSpinleRef' filesep 'kin_' num2str(kIdx,'%04d') '_P' num2str(poleIdx)];
        three=[imread([outputDirSlices1 filesep 'frame_nb' num2str(frameIdx,'%04d') '.png'])];
%     end
    writeVideo(video,three);
    %     fprintf('\b|\n');
end
close(video)
