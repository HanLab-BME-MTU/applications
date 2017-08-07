function projImages=project1D(MD,varargin)
% tracks are in original FoF and projected  in the manifold
% described by polygon.
% Poligon described the projected space (without dilated boundaries), It is composed of Tracks.
% in 1D, it is only composed of two tracks.
% dynPoligonISO define the ROI in the original, pixel ref and is only used for pixel masking (1D only).
% dynPoligonREF define the ROI in the final ref and is used for the actual displayed ROI location.
ip = inputParser;
ip.CaseSensitive = false;
ip.KeepUnmatched = true;
ip.addOptional('dynPoligonISO',[]);
ip.addOptional('dynPoligonREF',[]);
ip.addOptional('tracks',[]);
ip.addOptional('colormap',[]);
ip.addOptional('colorIndx',[]);
ip.addOptional('crop','manifold');
ip.addOptional('transType','affineOnePass');
ip.addOptional('processFrame',1:MD.nFrames_);
ip.addOptional('FoF',[]);
ip.addOptional('channelRender','grayRed');
ip.addOptional('intMinPrctil',[1 99.9]);
ip.addOptional('intMaxPrctil',[100 100]);
ip.addOptional('name',[]);
ip.addOptional('processSingleProj',[]);
ip.addOptional('fringeWidth',[]);
ip.addOptional('insetFringeWidth',20);
ip.parse(varargin{:});
p=ip.Results;
dynPoligonISO=p.dynPoligonISO;
dynPoligonREF=p.dynPoligonREF;
if(isempty(dynPoligonREF))
    dynPoligonREF=dynPoligonISO;
end
tracks=p.tracks;

showDebugGraphics=0;
insetFringeWidth=p.insetFringeWidth;
if(isempty(p.fringeWidth))
    fringeWidth=insetFringeWidth;
else
    fringeWidth=p.fringeWidth;
end
processFrame=p.processFrame;


%% Set normalization value
minIntensityNorm=[];
maxIntensityNorm=[];
for chIdx=1:length(MD.channels_)
    vol=MD.getChannel(chIdx).loadStack(1); 
    minIntensityNorm=[ minIntensityNorm prctile(vol(:),p.intMinPrctil(chIdx))];
    maxIntensityNorm=[ maxIntensityNorm prctile(vol(:),p.intMaxPrctil(chIdx))];
end

%% Define the static Rectangular cuboid that contains the pixel to be projected in the frame of reference.
%% Accordinly,  the coordinate of this cube are specified such as the origin of the frame of reference is the zero.

%% In the manifold crop case, the boundaries are given by the transform coordinate along the manifold polygon
%% in the full case, one have to estimate the maximum rectangle cuboid contained that can descibe the extremum coordinate of the original volume
if(strcmp(p.crop,'manifold')&&(~isempty(dynPoligonREF)))
    
    minX=MD.getDimensions('X')+1;
    minY=MD.getDimensions('Y')+1;
    minZ=(MD.pixelSizeZ_/MD.pixelSize_)*MD.getDimensions('Z')+1;
    maxX=0;
    maxY=0;
    maxZ=0;
    for iP=1:length(dynPoligonREF)
        minX=floor(min(min(dynPoligonREF(iP).x),minX));
        minY=floor(min(min(dynPoligonREF(iP).y),minY));
        minZ=floor(min(min(dynPoligonREF(iP).z),minZ));
        maxX=ceil(max(max(dynPoligonREF(iP).x),maxX));
        maxY=ceil(max(max(dynPoligonREF(iP).y),maxY));
        maxZ=ceil(max(max(dynPoligonREF(iP).z),maxZ));
    end
    maxXBorder=(maxX+fringeWidth);
    maxYBorder=(maxY+fringeWidth);
    maxZBorder=(maxZ+fringeWidth);
    minXBorder=(minX-fringeWidth);
    minYBorder=(minY-fringeWidth);
    minZBorder=(minZ-fringeWidth);
else
    if(~isempty(p.FoF))
        maxXBorder=MD.getDimensions('X')-p.FoF.origin(1,1);
        maxYBorder=MD.getDimensions('Y')-p.FoF.origin(1,2);
        maxZBorder=MD.getDimensions('Z')*(MD.pixelSizeZ_/MD.pixelSize_)-p.FoF.origin(1,3);
        minXBorder=1-p.FoF.origin(1,1);
        minYBorder=1-p.FoF.origin(1,2);
        minZBorder=1-p.FoF.origin(1,3);
    else
        maxXBorder=MD.getDimensions('X');
        maxYBorder=MD.getDimensions('Y');
        maxZBorder=MD.getDimensions('Z')*(MD.pixelSizeZ_/MD.pixelSize_);
        minXBorder=1;
        minYBorder=1;
        minZBorder=1;
    end
end

tform=affine3d();

if(~isempty(p.FoF))
    B=p.FoF.getBase(1);
    tform.T(1:3,1:3)=B;
end


%% create projection process saving independant projection location
outputDirSlices1=[MD.outputDirectory_ filesep '1DProjection' filesep p.name  ];
outputDirSingleProj=[MD.outputDirectory_ filesep '1DProjection' filesep p.name  ];
if(~isempty(p.processSingleProj))
    mkdirRobust([outputDirSingleProj filesep 'XY'])
    mkdirRobust([outputDirSingleProj filesep 'YZ'])
    mkdirRobust([outputDirSingleProj filesep 'XZ'])    
    p.processSingleProj.setOutFilePaths({[outputDirSingleProj filesep 'XY' filesep 'frame_nb%04d.png'], ...
        [outputDirSingleProj filesep 'YZ' filesep 'frame_nb%04d.png'], ...
        [outputDirSingleProj filesep 'XZ' filesep 'frame_nb%04d.png'], ...
        [outputDirSingleProj filesep 'limits.mat']});
    frameNb=MD.nFrames_;
    save([outputDirSingleProj filesep 'limits.mat'],'minXBorder', 'maxXBorder','minYBorder','maxYBorder','minZBorder','maxZBorder','frameNb');
end
outputDirDemo=[MD.outputDirectory_ filesep '1DProjection' filesep p.name filesep 'volDemo' ];

mkdirRobust(outputDirSlices1);
mkdirRobust(outputDirDemo);


% warp, crop, fuse and save each time point
parfor fIdx=processFrame
    disp(num2str(fIdx))
    
    % produce a ROI mask using the 1D polygon (segment defined by the extremities of the dynPoligonISO).
    % todo: N Channel (now 2).
    if(~isempty(dynPoligonISO))
         % Collect relative frameIdx
        pIndices=nan(1,length(dynPoligonISO));
        for polIdx=1:length(dynPoligonISO)
            F=dynPoligonISO(polIdx).f;
            pIdx=find(F==fIdx,1);
            if isempty(pIdx)
                if(fIdx>max(F))   pIdx=length(F);  else   pIdx=1; end;
            end
            pIndices(polIdx)=pIdx;
        end;

        %% Building mask in the 1D case
        nextPoint=length(dynPoligonISO);
        PCurrent=[dynPoligonISO(1).x(pIndices(1)) dynPoligonISO(1).y(pIndices(1)) dynPoligonISO(1).z(pIndices(1))];
        KCurrent=[dynPoligonISO(nextPoint).x(pIndices(nextPoint)) dynPoligonISO(nextPoint).y(pIndices(nextPoint)) dynPoligonISO(nextPoint).z(pIndices(nextPoint))];

        % Building mask for both channel on the whole volume
        % NOTE: in order to apply fringe isotropically, we need the mask to
        % be isotropized briefly.
        mask=zeros(size(vol,1),size(vol,2),ceil(size(vol,3)*MD.pixelSizeZ_/MD.pixelSize_));
        sampling=100;
        xSeg=round(linspace(PCurrent(1),KCurrent(1),sampling));
        ySeg=round(linspace(PCurrent(2),KCurrent(2),sampling));
        zSeg=round(linspace(PCurrent(3),KCurrent(3),sampling));
        indx=sub2ind(size(mask),ySeg,xSeg,zSeg);
        
        mask(indx)=1;
        
        %mask=imdilate(mask,IMSphere);  %ones(cubeHalfWidth,cubeHalfWidth,round(cubeHalfWidth*MD.pixelSize_/MD.pixelSizeZ_)));
        %% If no transform are needed, now to save on bwdist.

        distMap=mask;

%         mask=resize(distMap,[1 1 MD.pixelSizeZ_/MD.pixelSize_]);
        distMap=bwdist(distMap);
%         distMap=resize(distMap,[1 1 MD.pixelSize_/MD.pixelSizeZ_]);
        mask(distMap<insetFringeWidth)=1;
        [y x z]=...
            ndgrid( linspace(1,size(mask,1),size(vol,1)),...
                    linspace(1,size(mask,2),size(vol,2)),...
                    linspace(1,size(mask,3),size(vol,3)));
        mask=interp3(mask,x,y,z);
    end
    mips=cell(3,length(MD.channels_));
   
    for chIdx=1:length(MD.channels_)
        vol=MD.getChannel(chIdx).loadStack(fIdx);
        if(~isempty(dynPoligonISO))
            if(isempty(p.FoF))&&(strcmp(p.crop,'manifold'))
                aminXBorder=max(1,minXBorder);
                aminYBorder=max(1,minYBorder);
                aminZBorder=max(1,minZBorder);
                amaxXBorder=min(size(mask,2),maxXBorder);
                amaxYBorder=min(size(mask,1),maxYBorder);
                amaxZBorder=min(size(mask,3),maxZBorder);
                
                XCropMask=false(1,size(mask,2));
                XCropMask(aminXBorder:amaxXBorder)=true;
                YCropMask=false(1,size(mask,1));
                YCropMask(aminYBorder:amaxYBorder)=true;
                ZCropMask=false(1,size(mask,3));
                ZCropMask(ceil(aminZBorder:amaxZBorder))=true;
                ZCropVol=false(1,size(vol,3));
                ZCropVol(ceil((aminZBorder:amaxZBorder)*MD.pixelSize_/MD.pixelSizeZ_))=true;
                mask(:,:,~ZCropMask)=[];
                mask(~YCropMask,:,:)=[];
                mask(:,~XCropMask,:)=[];
                
                vol(:,:,~ZCropVol)=[];
                vol(~YCropMask,:,:)=[];
                vol(:,~XCropMask,:)=[];
                
            end
            
            maskedVol=vol;
            maskedVol(~mask)=0;
        end
        
        % If needed the map must rotated before cropped (efficiency)
        % Rotation depends on the FrameOfRef associated to the tracks the compose the dynanimc polygon
        % Cropping area according to the polygon OVER TIME plus added vizualiation margin
        % Rotation will use imwarp
        % Can we use imwar for cropping too ?
        
        %% if a FoF is specified, warp and crop data according to the
        tform=affine3d();
        warpedVol=vol;
        if(~isempty(dynPoligonISO))
            warpedMaskedVol=maskedVol;
        else
            warpedMaskedVol=zeros(size(warpedVol));
        end
        if(~isempty(p.FoF))
            B=p.FoF.getBase(fIdx);
            tform.T(4,[1 2 3])=(-p.FoF.getOrigAtFrame(fIdx)+p.FoF.origin(1,:))*B;
            tform.T(1:3,1:3)=B;
            %
            tformTransOnly=affine3d();
            tformTransOnly.T(4,[1 2 3])=(-p.FoF.getOrigAtFrame(fIdx));
            
            %
            tformRelTransOnly=affine3d();
            tformRelTransOnly.T(4,[1 2 3])=(-p.FoF.origin(1,:)+p.FoF.getOrigAtFrame(fIdx));
            
            tformRotOnly=affine3d();
            B=p.FoF.getBase(fIdx);
            tformRotOnly.T(1:3,1:3)=B;
            
            tformRotOnlyInit=affine3d();
            B=p.FoF.getBase(1);
            tformRotOnlyInit.T(1:3,1:3)=B;
            
            orig=p.FoF.getOrigAtFrame(fIdx);
            
            inputRef=imref3d([ MD.getDimensions('Y') MD.getDimensions('X') MD.getDimensions('Z')], ...
                [1 MD.getDimensions('X')],[1 MD.getDimensions('Y')],[1 MD.getDimensions('Z')*MD.pixelSizeZ_/MD.pixelSize_]);
            
            
            switch p.transType
                case 'affineOnePass'
                    maxXBorderFull=MD.getDimensions('X');
                    maxYBorderFull=MD.getDimensions('Y');
                    maxZBorderFull=MD.getDimensions('Z')*(MD.pixelSizeZ_/MD.pixelSize_);
                    minXBorderFull=1;
                    minYBorderFull=1;
                    minZBorderFull=1;
                    
                    minXBorderCurr=minXBorderFull - orig(1); maxXBorderCurr=maxXBorderFull -  orig(1);
                    minYBorderCurr=minYBorderFull - orig(2); maxYBorderCurr=maxYBorderFull  - orig(2);
                    minZBorderCurr=minZBorderFull - orig(3); maxZBorderCurr=maxZBorderFull -  orig(3);
                    
                    inputRef=imref3d([ MD.getDimensions('Y') MD.getDimensions('X') MD.getDimensions('Z')], ...
                        [minXBorderCurr maxXBorderCurr], [minYBorderCurr maxYBorderCurr], [minZBorderCurr maxZBorderCurr]);
                    
                    minXBorderCurr=minXBorder ;%+ orig(1) - p.FoF.origin(1,1);
                    maxXBorderCurr=maxXBorder ;%+ orig(1) - p.FoF.origin(1,1);
                    minYBorderCurr=minYBorder ;%+ orig(2) - p.FoF.origin(1,2);
                    maxYBorderCurr=maxYBorder ;%+ orig(2) - p.FoF.origin(1,2);
                    minZBorderCurr=minZBorder ;%+ orig(3) - p.FoF.origin(1,3);
                    maxZBorderCurr=maxZBorder ;%+ orig(3) - p.FoF.origin(1,3);
                    
                    rotOutputRef=imref3d([    ceil(maxYBorderCurr-minYBorderCurr) ...
                        ceil(maxXBorderCurr-minXBorderCurr) ...
                        ceil(maxZBorderCurr-minZBorderCurr) ], ...
                        [minXBorderCurr maxXBorderCurr], [minYBorderCurr maxYBorderCurr], [minZBorderCurr maxZBorderCurr]);
                    
                    
                    warpedVol=imwarp(vol,inputRef,tformRotOnly,'OutputView',rotOutputRef);
                    
                    if(~isempty(dynPoligonISO))
                        warpedMaskedVol=imwarp(maskedVol,inputRef,tformRotOnly,'OutputView',rotOutputRef);
                    else
                        warpedMaskedVol=zeros(size(warpedVol));
                    end
                    
                case 'translation'
                    disp(num2str(fIdx))
                    minXBorderCurr=minXBorder ;%+ orig(1) - p.FoF.origin(1,1);
                    maxXBorderCurr=maxXBorder ;%+ orig(1) - p.FoF.origin(1,1);
                    minYBorderCurr=minYBorder ;%+ orig(2) - p.FoF.origin(1,2);
                    maxYBorderCurr=maxYBorder ;%+ orig(2) - p.FoF.origin(1,2);
                    minZBorderCurr=minZBorder ;%+ orig(3) - p.FoF.origin(1,3);
                    maxZBorderCurr=maxZBorder ;%+ orig(3) - p.FoF.origin(1,3);
                    
                    %             [xLimitsOut,yLimitsOut,zLimitsOut] = outputLimits(tformTransOnly,[minXBorder maxXBorder], [minYBorder maxYBorder], [minZBorder maxZBorder]);
                    %             minXBorderCurr=xLimitsOut(1); maxXBorderCurr=xLimitsOut(2);
                    %             minYBorderCurr=yLimitsOut(1); maxYBorderCurr=yLimitsOut(2);
                    %             minZBorderCurr=zLimitsOut(1); maxZBorderCurr=zLimitsOut(2);
                    %
                    outputRef=imref3d([ ceil(maxYBorderCurr-minYBorderCurr) ...
                        ceil(maxXBorderCurr-minXBorderCurr) ...
                        ceil(maxZBorderCurr-minZBorderCurr) ], ...
                        [minXBorderCurr maxXBorderCurr], [minYBorderCurr maxYBorderCurr], [minZBorderCurr maxZBorderCurr]);
                    
                    warpedVol=imwarp(vol,inputRef,tformTransOnly,'OutputView',outputRef);
                    warpedMaskedVol=imwarp(maskedVol,inputRef,tformTransOnly,'OutputView',outputRef);
                    
                case 'transCrop'
                    %% to be updated
                    [xLimitsOut,yLimitsOut,zLimitsOut] = outputLimits(tformRelTransOnly,[minXBorder maxXBorder], [minYBorder maxYBorder], [minZBorder maxZBorder]);
                    minXBorderCurr=xLimitsOut(1); maxXBorderCurr=xLimitsOut(2);
                    minYBorderCurr=yLimitsOut(1); maxYBorderCurr=yLimitsOut(2);
                    minZBorderCurr=zLimitsOut(1); maxZBorderCurr=zLimitsOut(2);
                    
                    maskcrop=maskedVol;
                    nullMaskXY=(squeeze(any(maskcrop,3)));
                    YNull=~(squeeze(any(any(mask,3),2)));
                    XNull=~(squeeze(any(any(mask,3),1)));
                    ZNull=~(squeeze(any(any(mask,1),2)));
                    
                    YNull= zeros(1,size(maskcrop,1));
                    YNull(1:minYBorderCurr)=1;
                    YNull(maxYBorderCurr:end)=1;
                    YNull=logical(YNull);
                    
                    XNull= zeros(1,size(maskcrop,2));
                    XNull(1:minXBorderCurr)=1;
                    XNull(maxXBorderCurr:end)=1;
                    XNull=logical(XNull);
                    
                    ZNull= zeros(1,size(maskcrop,3));
                    ZNull(1:ceil(minZBorderCurr*MD.pixelSize_/MD.pixelSizeZ_))=1;
                    ZNull(ceil(maxZBorderCurr*MD.pixelSize_/MD.pixelSizeZ_):end)=1;
                    ZNull=logical(ZNull);
                    
                    maskcrop(:,:,ZNull)=[];
                    maskcrop(YNull,:,:)=[];
                    maskcrop(:,XNull,:)=[];
                    
                    warpedMaskedVol=maskcrop;
                    
                    warpedVol=vol;
                    warpedVol(:,:,ZNull)=[];
                    warpedVol(YNull,:,:)=[];
                    warpedVol(:,XNull,:)=[];
                    
                otherwise
                    error('unknown trans type');
            end
        end
        
        
        %% Create MIPS for each channel, fuse mask and full volume
        switch p.transType
            case 'none'
            case 'transCrop'
                ZRatio=MD.pixelSizeZ_/MD.pixelSize_;
            otherwise
                ZRatio=1;
        end;
        if(isempty(p.FoF))
            ZRatio=MD.pixelSizeZ_/MD.pixelSize_;
        end
        % Create MIP of ROI and context
        [fullmaxXY,fullmaxZY,fullmaxZX,~]=computeMIPs(warpedVol,ZRatio, ...
            minIntensityNorm(chIdx),maxIntensityNorm(chIdx));
        [maxXY,maxZY,maxZX,~]=computeMIPs(warpedMaskedVol,ZRatio, ...
            minIntensityNorm(chIdx),maxIntensityNorm(chIdx));
        
        % Fuse ROI and context
        maxXY(maxXY==0)=fullmaxXY(maxXY==0);
        maxZY(maxZY==0)=fullmaxZY(maxZY==0);
        maxZX(maxZX==0)=fullmaxZX(maxZX==0);
        
        
        %% Resize and fuse channel MIPS
        maxMIPSize=400;
        [sX,sY,sZ]=size(warpedMaskedVol);
        if(strcmp(p.transType,'none'))
            sZ=sZ*MD.pixelSizeZ_/MD.pixelSize_;
        end
        resizeScale=maxMIPSize/max([sX,sY,sZ]);
        
        XYMax=imresize(maxXY,resizeScale);
        ZYMax=imresize(maxZY,resizeScale);
        ZXMax=imresize(maxZX,resizeScale);
        mips{1,chIdx}=XYMax;
        mips{2,chIdx}=ZYMax;
        mips{3,chIdx}=ZXMax;
    end
    
    % fuse volume if 2 channels
    if(length(MD.channels_)==2)
        XYProj=renderChannel(mips{1,1},mips{1,2},p.channelRender);
        ZYProj=renderChannel(mips{2,1},mips{2,2},p.channelRender);
        ZXProj=renderChannel(mips{3,1},mips{3,2},p.channelRender);
    else
        XYProj=repmat(mips{1,1},1,1,3);
        ZYProj=repmat(mips{2,1},1,1,3);
        ZXProj=repmat(mips{3,1},1,1,3); 
    end
    
    %% select tracks that starts or end in the manifold
    inMaskTrack=zeros(1,length(tracks));
    for tIdx=1:length(tracks)
        t=tracks(tIdx);
        pIdx=find(t.f==fIdx,1);
        if(~isempty(pIdx))
            Vq = interp3(   linspace(minXBorder,maxXBorder,size(warpedMaskedVol,1)), ...
                linspace(minYBorder,maxYBorder,size(warpedMaskedVol,2)), ...
                linspace(minZBorder,maxZBorder,size(warpedMaskedVol,3)), ...
                warpedMaskedVol,t.x(pIdx), t.y(pIdx), t.z(pIdx));
            inMaskTrack(tIdx)=Vq>0;
        end
    end
    tracksInMask=tracks(logical(inMaskTrack));

    [tracksXY,tracksZY,tracksZX]=overlayProjTracks(XYProj,ZYProj,ZXProj, ...
      [minXBorder maxXBorder],[minYBorder maxYBorder],[minZBorder maxZBorder], ...
      fIdx,tracksInMask,p.colormap,p.colorIndx(logical(inMaskTrack)));

    %% write images
    if(~isempty(p.processSingleProj))
      imwrite(tracksXY,[outputDirSingleProj filesep 'XY' filesep  'frame_nb' num2str(fIdx,'%04d') '.png']);
      imwrite(tracksZY,[outputDirSingleProj filesep 'YZ' filesep  'frame_nb' num2str(fIdx,'%04d') '.png']);
      imwrite(tracksZX,[outputDirSingleProj filesep 'XZ' filesep  'frame_nb' num2str(fIdx,'%04d') '.png']);
    end

    %% Use Z to index image line (going up)
    tracksZY=permute(tracksZY,[2 1 3]);
    tracksZX=permute(tracksZX,[2 1 3]);
    three=projMontage(tracksXY,tracksZX,tracksZY);
    imwrite(three,[outputDirSlices1 filesep 'frame_nb' num2str(fIdx,'%04d') '.png']);
end

video = VideoWriter([outputDirSlices1  '.avi']);
video.FrameRate = 5;  % Default 30
video.Quality = 100;    % Default 75

open(video)
for frameIdx=processFrame
    % save the maximum intensity projections
%     for poleIdx=poleIndices
%         outputDirSlices1=[MD.outputDirectory_ filesep 'maskSliceTracksSpinleRef' filesep 'kin_' num2str(kIdx,'%04d') '_P' num2str(poleIdx)];
        three=[imread([outputDirSlices1 filesep 'frame_nb' num2str(frameIdx,'%04d') '.png'])];
%     end
    writeVideo(video,three);
    %     fprintf('\b|\n');
end
close(video)


function RGBVol=renderChannel(ch1,ch2,type)
    switch(type)
        case 'grayRed'
          RGBVol=grayRedRender(ch1,ch2);
        case 'greenRed'
            RGBVol=greenRedRender(ch1,ch2);
        otherwise
            error('unknown channel renderer');
    end

function RGBVol=greenRedRender(greenCh,redCh)
    RGBVol=repmat(greenCh,1,1,3);
    %RGBThree(:,:,1)=max(rMaxXY,rmaxXYKin);
    RGBVol(:,:,1)=redCh;
    RGBVol(:,:,2)=greenCh;
    RGBVol(:,:,3)=0;


function RGBVol=grayRedRender(grayCh,redCh)
    RGBVol=repmat(grayCh,1,1,3);
    RGBVol(:,:,1)=max(grayCh,redCh);
    
    

