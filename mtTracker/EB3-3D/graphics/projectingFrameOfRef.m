function projImages=project1D(MD,frameOfRef,dynPoligon,tracks,tracksColors)
% FrameOfRef is the frameOfRef in which resides the tracks that are to be projected in the manifold
% described by polygon.
% Poligon described the projected space (without dilated boundaries), It is composed of Tracks.
% in 1D, it is only composed of two tracks.

%% Test on a single kinetochore bundle

% outputDirBundle=[MD.outputDirectory_ filesep 'Kin' filesep 'bundles'];
% tmp=load([outputDirBundle filesep 'kin-MT-bundle.mat'],'kinTracks');
% kinTracksBundle=tmp.kinTracks;
%%

showDebugGraphics=0;
cubeHalfWidth=20;

system(['mkdir -p ' outputDirSlices1]);
outputDirSlices1=[MD.outputDirectory_ filesep '1DProjection' filesep 'kin_' num2str(kIdx,'%04d') '_P' num2str(poleIdx)];
for fIdx=1:MD.nFrames_
  vol=MD.getChannel(1).loadStack(fIdx);
  kinvol=MD.getChannel(2).loadStack(fIdx);
  pIndices=arrayfun(@(dp) find(dp.f==fIdx),dynPoligon);
  if(~isempty(kinTIdx))
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
    %stackWrite(uint16(maskedVol),[outputDirDemo filesep 'volDemo'  num2str(fIdx,'%04d') '.tif']) % Names the file stitched001..., in /stitched/

    maskedKin=kinvol;
    maskedKin(~mask)=0;
    stackWrite(uint16(maskedKin),[outputDirDemoKin filesep 'volDemoKin'  num2str(fIdx,'%04d') '.tif']) % Names the file stitched001..., in /stitched/

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


    resizeScale=400/size(maxXY,1);
    rMaxXY=imresize(maxXY,resizeScale);
    rmaxXYKin=imresize(maxXYKin,resizeScale);

    RGBThree=repmat(rMaxXY,1,1,3);
    % channel fusion
    RGBThree(:,:,1)=max(rMaxXY,rmaxXYKin);

    myColormap=[[0 0 255];
    [0 255 00]];

    tracksColors=uint8(myColormap);
    tracksXY=trackBinaryOverlay(RGBThree,[find(~XNull,1) find(~XNull,1,'last')],[find(~YNull,1) find(~YNull,1,'last')],capturedEB3OrigRef,fIdx,(kinTracksBundle(kIdx).fiber>0)+1,tracksColors);
    %%
    rMax=imresize(maxZY,[mipSize(1) mipSize(2)*MD.pixelSize_/MD.pixelSizeZ_]);
    rmaxKin=imresize(maxZYKin,[mipSize(1) mipSize(2)*MD.pixelSize_/MD.pixelSizeZ_]);

    RGBThree=repmat(rMax,1,1,3);
    RGBThree(:,:,1)=max(rMax,rmaxKin);

    capturedEB3ZY=capturedEB3OrigRef.copy();
    for ebIdx=1:length(capturedEB3ZY)
      capturedEB3ZY(ebIdx).x=capturedEB3OrigRef(ebIdx).z*MD.pixelSize_/MD.pixelSizeZ_;
    end
    tracksColors=uint8(myColormap);
    tracksZY=trackBinaryOverlay(RGBThree,[find(~ZNull,1) find(~ZNull,1,'last')],[find(~YNull,1) find(~YNull,1,'last')],capturedEB3ZY,fIdx,(kinTracksBundle(kIdx).fiber>0)+1,tracksColors);

    tracksZY=permute(tracksZY,[2 1 3]);
    stripeSize=4;
    threeTop = [tracksXY, 0*ones(size(tracksXY,1), stripeSize,3), 0*ones(size(tracksXY,1), size(tracksZY,2),3)];

    %%
    rMax=imresize(maxZX,[mipSize(1) mipSize(2)*MD.pixelSize_/MD.pixelSizeZ_]);
    rmaxKin=imresize(maxZXKin,[mipSize(1) mipSize(2)*MD.pixelSize_/MD.pixelSizeZ_]);

    RGBThree=repmat(rMax,1,1,3);
    RGBThree(:,:,1)=max(rMax,rmaxKin);

    capturedEB3ZX=capturedEB3ZY.copy();
    for ebIdx=1:length(capturedEB3ZX)
      capturedEB3ZX(ebIdx).y=capturedEB3OrigRef(ebIdx).x;
    end
    tracksColors=uint8(myColormap);
    tracksZX=trackBinaryOverlay(RGBThree,[find(~ZNull,1) find(~ZNull,1,'last')],[find(~XNull,1) find(~XNull,1,'last')],capturedEB3ZX,fIdx,(kinTracksBundle(kIdx).fiber>0)+1,tracksColors);

    tracksZX=permute(tracksZX,[2 1 3]);
    threeBottom = [tracksZX, 0*ones(size(tracksZX,1),+stripeSize,3),tracksZY];
    three = [threeTop; ones(stripeSize, size(tracksXY,2)+size(tracksZY,2)+stripeSize,3); threeBottom];
    %%
    imwrite(three,[outputDirSlices1 filesep 'frame_nb' num2str(fIdx,'%04d') '.png']);
  end
end
end
video = VideoWriter([outputDirSlices1  '.avi']);
video.FrameRate = 5;  % Default 30
video.Quality = 100;    % Default 75

open(video)
for frameIdx=1:150
    % save the maximum intensity projections
    three=[];
    for poleIdx=poleIndices
        outputDirSlices1=[MD.outputDirectory_ filesep 'maskSliceTracksSpinleRef' filesep 'kin_' num2str(kIdx,'%04d') '_P' num2str(poleIdx)];
        three=[three imread([outputDirSlices1 filesep 'frame_nb' num2str(frameIdx,'%04d') '.png'])];
    end
    writeVideo(video,three);
%     fprintf('\b|\n');
end
close(video)
