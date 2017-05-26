function projImages=project1D(MD,dynPoligonISO,varargin)
% tracks are in original FoF and projected  in the manifold
% described by polygon.
% Poligon described the projected space (without dilated boundaries), It is composed of Tracks.
% in 1D, it is only composed of two tracks.
% dynPoligonISO define the ROI in the original, pixel ref and is only used for pixel masking (1D only).
% dynPoligonREF define the ROI in the final ref and is used for the actual displayed ROI location.
ip = inputParser;
ip.CaseSensitive = false;
ip.KeepUnmatched = true;
ip.addOptional('dynPoligonREF',dynPoligonISO);
ip.addOptional('tracks',[]);
ip.addOptional('colormap',[]);
ip.addOptional('colorIndx',[]);
ip.addOptional('crop','manifold');
ip.addOptional('transType','affineOnePass');
ip.addOptional('FoF',[]);
ip.addOptional('channelRender','grayRed');
ip.addOptional('intMinPrctil',[1 99.9]);
ip.addOptional('intMaxPrctil',[100 100]);
ip.addOptional('name',[]);
ip.addOptional('processSingleProj',[]);
ip.addOptional('fringeWidth',20);
ip.parse(varargin{:});
p=ip.Results;

tracks=p.tracks;

showDebugGraphics=0;
fringeWidth=p.fringeWidth;

%% Define the static Rectangular cuboid that contains the pixel to be projected in the frame of reference.
%% Accordinly,  the coordinate of this cube are specified such as the origin of the frame of reference is the zero.

%% In the manifold crop case, the boundaries are given by the transform coordinate along the manifold polygon
%% in the full case, one have to estimate the maximum rectangle cuboid contained that can descibe the extremum coordinate of the original volume
orig=p.FoF.origin;
dynPoligonREF=p.dynPoligonREF;
if(strcmp(p.crop,'manifold'))
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
    % minZ=floor(minZ*MD.pixelSize_/MD.pixelSizeZ_);
    % maxZ=ceil(maxZ*MD.pixelSize_/MD.pixelSizeZ_);
%
%     maxXBorder=min(maxX+cubeHalfWidth,MD.getDimensions('X'));
%     maxYBorder=min(maxY+cubeHalfWidth,MD.getDimensions('Y'));
%     maxZBorder=min(maxZ+cubeHalfWidth,(MD.pixelSizeZ_/MD.pixelSize_)*MD.getDimensions('Z'));
%     minXBorder=max(minX-cubeHalfWidth,1);
%     minYBorder=max(minY-cubeHalfWidth,1);
%     minZBorder=max(minZ-cubeHalfWidth,1);

    maxXBorder=(maxX+fringeWidth);
    maxYBorder=(maxY+fringeWidth);
    maxZBorder=(maxZ+fringeWidth);
    minXBorder=(minX-fringeWidth);
    minYBorder=(minY-fringeWidth);
    minZBorder=(minZ-fringeWidth);
else
    maxXBorder=MD.getDimensions('X')-p.FoF.origin(1,1);
    maxYBorder=MD.getDimensions('Y')-p.FoF.origin(1,2);
    maxZBorder=MD.getDimensions('Z')*(MD.pixelSizeZ_/MD.pixelSize_)-p.FoF.origin(1,3);
    minXBorder=1-p.FoF.origin(1,1);
    minYBorder=1-p.FoF.origin(1,2);
    minZBorder=1-p.FoF.origin(1,3);
end

tform=affine3d();

if(~isempty(p.FoF))
    B=p.FoF.getBase(1);
    tform.T(1:3,1:3)=B;
end


%% testing imwarp to crop the image
outputDirSlices1=[MD.outputDirectory_ filesep '1DProjection' filesep p.name  ];
outputDirSingleProj=[MD.outputDirectory_ filesep '1DProjection' filesep p.name  ];
if(~isempty(p.processSingleProj))
    mkdirRobust([outputDirSingleProj filesep 'XY'])
    mkdirRobust([outputDirSingleProj filesep 'XZ'])
    mkdirRobust([outputDirSingleProj filesep 'YZ'])
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
parfor fIdx=1:MD.nFrames_
    disp(num2str(fIdx))

    vol=MD.getChannel(1).loadStack(fIdx);
    kinvol=MD.getChannel(2).loadStack(fIdx);
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
        % NOTE: Could masking use imwarp for speed ?
        mask=zeros(size(vol,1),size(vol,2),ceil(size(vol,3)*MD.pixelSizeZ_/MD.pixelSize_));
        sampling=100;
        xSeg=round(linspace(PCurrent(1),KCurrent(1),sampling));
        ySeg=round(linspace(PCurrent(2),KCurrent(2),sampling));
        zSeg=round(linspace(PCurrent(3),KCurrent(3),sampling));
        indx=sub2ind(size(mask),ySeg,xSeg,zSeg);
        
        mask(indx)=1;
        %mask=imdilate(mask,IMSphere);  %ones(cubeHalfWidth,cubeHalfWidth,round(cubeHalfWidth*MD.pixelSize_/MD.pixelSizeZ_)));
        distMap=mask;
%         mask=resize(distMap,[1 1 MD.pixelSizeZ_/MD.pixelSize_]);
        distMap=bwdist(distMap);
%         distMap=resize(distMap,[1 1 MD.pixelSize_/MD.pixelSizeZ_]);
        mask(distMap<fringeWidth)=1;
        [y x z]=...
            ndgrid( linspace(1,size(mask,1),size(vol,1)),...
                    linspace(1,size(mask,2),size(vol,2)),...
                    linspace(1,size(mask,3),size(vol,3)));
        mask=interp3(mask,x,y,z);
        maskedVol=vol;
        maskedVol(~mask)=0;
        %     outputDir=[outputDirDemo filesep  'mask'];mkdir(outputDir);
        %     imwrite(uint8(255*mat2gray(max(maskedVol,[],3))),[outputDir filesep  'mask_' num2str(fIdx,'%04d') '.png']) % Names the file stitched001..., in /stitched/

        maskedKin=kinvol;
        maskedKin(~mask)=0;
        %     outputDir=[outputDirDemo filesep  'maskKin_'];mkdir(outputDir);
        %     imwrite(uint8(255*mat2gray(max(maskedKin,[],3))),[outputDir filesep  'maskKin_' num2str(fIdx,'%04d') '.png']) % Names the file stitched001..., in /stitched/
    end

    % If needed the map must rotated before cropped (efficiency)
    % Rotation depends on the FrameOfRef associated to the tracks the compose the dynanimc polygon
    % Cropping area according to the polygon OVER TIME plus added vizualiation margin
    % Rotation will use imwarp
    % Can we use imwar for cropping too ?

    tform=affine3d();

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
    end
    orig=p.FoF.getOrigAtFrame(fIdx);

    inputRef=imref3d([ MD.getDimensions('Y') MD.getDimensions('X') MD.getDimensions('Z')], ...
        [1 MD.getDimensions('X')],[1 MD.getDimensions('Y')],[1 MD.getDimensions('Z')*MD.pixelSizeZ_/MD.pixelSize_]);


    switch p.transType
        case 'affine'
            maxXBorderFull=MD.getDimensions('X');
            maxYBorderFull=MD.getDimensions('Y');
            maxZBorderFull=MD.getDimensions('Z')*(MD.pixelSizeZ_/MD.pixelSize_);
            minXBorderFull=1;
            minYBorderFull=1;
            minZBorderFull=1;

            border=300;
            minXBorderCurr=minXBorderFull - orig(1)-border; maxXBorderCurr=maxXBorderFull - orig(1)+border;
            minYBorderCurr=minYBorderFull - orig(2)-border; maxYBorderCurr=maxYBorderFull - orig(2)+border;
            minZBorderCurr=minZBorderFull - orig(3)-border; maxZBorderCurr=maxZBorderFull - orig(3)+border;

            outputRef=imref3d([ ceil(maxYBorderFull-minYBorderFull) ...
                                ceil(maxXBorderFull-minXBorderFull) ...
                                ceil(maxZBorderFull-minZBorderFull)], ...
                                [minXBorderCurr maxXBorderCurr], ...
                                [minYBorderCurr maxYBorderCurr], ...
                                [minZBorderCurr maxZBorderCurr]);

            transVol=imwarp(vol,inputRef,tformTransOnly,'OutputView',outputRef);
            transMaskedVol=imwarp(maskedVol,inputRef,tformTransOnly,'OutputView',outputRef);

            transKinVol=imwarp(kinvol,inputRef,tformTransOnly,'OutputView',outputRef);
            transMaskedKinVol=imwarp(maskedKin,inputRef,tformTransOnly,'OutputView',outputRef);

            disp(num2str(fIdx))
            minXBorderCurr=minXBorder ;%+ orig(1) - p.FoF.origin(1,1);
            maxXBorderCurr=maxXBorder ;%+ orig(1) - p.FoF.origin(1,1);
            minYBorderCurr=minYBorder ;%+ orig(2) - p.FoF.origin(1,2);
            maxYBorderCurr=maxYBorder ;%+ orig(2) - p.FoF.origin(1,2);
            minZBorderCurr=minZBorder ;%+ orig(3) - p.FoF.origin(1,3);
            maxZBorderCurr=maxZBorder ;%+ orig(3) - p.FoF.origin(1,3);

%             [xLimitsOut,yLimitsOut,zLimitsOut] = outputLimits(tformRotOnlyInit,[minXBorderCurr maxXBorderCurr], [minYBorderCurr maxYBorderCurr], [minZBorderCurr maxZBorderCurr]);
%             minXBorderCurr=xLimitsOut(1); maxXBorderCurr=xLimitsOut(2);
%             minYBorderCurr=yLimitsOut(1); maxYBorderCurr=yLimitsOut(2);
%             minZBorderCurr=zLimitsOut(1); maxZBorderCurr=zLimitsOut(2);

            rotOutputRef=imref3d([ ceil(maxYBorderCurr-minYBorderCurr) ...
                ceil(maxXBorderCurr-minXBorderCurr) ...
                ceil(maxZBorderCurr-minZBorderCurr) ], ...
                [minXBorderCurr maxXBorderCurr], [minYBorderCurr maxYBorderCurr], [minZBorderCurr maxZBorderCurr]);

            warpedVol=imwarp(transVol,outputRef,tformRotOnly,'OutputView',rotOutputRef);
            warpedMaskedVol=imwarp(transMaskedVol,outputRef,tformRotOnly,'OutputView',rotOutputRef);

            warpedKinVol=imwarp(transKinVol,outputRef,tformRotOnly,'OutputView',rotOutputRef);
            warpedMaskedKinVol=imwarp(transMaskedKinVol,outputRef,tformRotOnly,'OutputView',rotOutputRef);

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

            % border=0;
            % minXBorderCurr=minXBorder + orig(1) - p.FoF.origin(1,1)-border;
            % maxXBorderCurr=maxXBorder + orig(1) -  p.FoF.origin(1,1)-border;
            % minYBorderCurr=minYBorder + orig(2) - p.FoF.origin(1,2)-border;
            % maxYBorderCurr=maxYBorder + orig(2) - p.FoF.origin(1,2)-border;
            % minZBorderCurr=minZBorder + orig(3) - p.FoF.origin(1,3)-border;
            % maxZBorderCurr=maxZBorder + orig(3) -  p.FoF.origin(1,3)-border;

            %% if full
            % [xLimitsOut,yLimitsOut,zLimitsOut] = outputLimits(tformRotOnlyInit,[minXBorderCurr maxXBorderCurr], [minYBorderCurr maxYBorderCurr], [minZBorderCurr maxZBorderCurr]);
            % minXBorderCurr=xLimitsOut(1); maxXBorderCurr=xLimitsOut(2);
            % minYBorderCurr=yLimitsOut(1); maxYBorderCurr=yLimitsOut(2);
            % minZBorderCurr=zLimitsOut(1); maxZBorderCurr=zLimitsOut(2);

            rotOutputRef=imref3d([    ceil(maxYBorderCurr-minYBorderCurr) ...
                                      ceil(maxXBorderCurr-minXBorderCurr) ...
                                      ceil(maxZBorderCurr-minZBorderCurr) ], ...
                [minXBorderCurr maxXBorderCurr], [minYBorderCurr maxYBorderCurr], [minZBorderCurr maxZBorderCurr]);


            warpedVol=imwarp(vol,inputRef,tformRotOnly,'OutputView',rotOutputRef);
            warpedKinVol=imwarp(kinvol,inputRef,tformRotOnly,'OutputView',rotOutputRef);

            if(~isempty(dynPoligonISO))
                warpedMaskedVol=imwarp(maskedVol,inputRef,tformRotOnly,'OutputView',rotOutputRef);
                warpedMaskedKinVol=imwarp(maskedKin,inputRef,tformRotOnly,'OutputView',rotOutputRef);
            else
                warpedMaskedVol=zeros(size(warpedVol));
                warpedMaskedKinVol=zeros(size(warpedKinVol));
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

                warpedKinVol=imwarp(kinvol,inputRef,tformTransOnly,'OutputView',outputRef);
                warpedMaskedKinVol=imwarp(maskedKin,inputRef,tformTransOnly,'OutputView',outputRef);

        case 'none'
            maskcrop=maskedVol;
            nullMaskXY=(squeeze(any(maskcrop,3)));
            YNull=~(squeeze(any(any(mask,3),2)));
            XNull=~(squeeze(any(any(mask,3),1)));
            ZNull=~(squeeze(any(any(mask,1),2)));

            minYBorderCurr=minYBorder + orig(2); % canonical ref projection
            YNull= zeros(1,size(maskcrop,1));
            YNull(1:minYBorderCurr)=1;
            YNull(minYBorderCurr:end)=1;
            YNull=logical(YNull);

            minXBorderCurr=minXBorder + orig(1); % canonical ref projection
            XNull= zeros(1,size(maskcrop,2));
            XNull(1:minXBorderCurr)=1;
            XNull(minXBorderCurr:end)=1;
            XNull=logical(XNull);

            minZBorderCurr=minZBorder + orig(3); % canonical ref projection
            ZNull= zeros(1,size(maskcrop,3));
            ZNull(1:ceil(minZBorderCurr*MD.pixelSize_/MD.pixelSizeZ_))=1;
            ZNull(ceil(minZBorderCurr*MD.pixelSize_/MD.pixelSizeZ_):end)=1;
            ZNull=logical(ZNull);

            maskcrop(:,:,ZNull)=[];
            maskcrop(YNull,:,:)=[];
            maskcrop(:,XNull,:)=[];
            maskedKin(:,:,ZNull)=[];
            maskedKin(YNull,:,:)=[];
            maskedKin(:,XNull,:)=[];

            warpedMaskedVol=maskcrop;
            warpedMaskedKinVol=maskedKin;

            warpedVol=vol;
            warpedKinVol=kinvol;
            warpedVol(:,:,ZNull)=[];
            warpedVol(YNull,:,:)=[];
            warpedVol(:,XNull,:)=[];
            warpedKinVol(:,:,ZNull)=[];
            warpedKinVol(YNull,:,:)=[];
            warpedKinVol(:,XNull,:)=[];

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
            maskedKin(:,:,ZNull)=[];
            maskedKin(YNull,:,:)=[];
            maskedKin(:,XNull,:)=[];

            warpedMaskedVol=maskcrop;
            warpedMaskedKinVol=maskedKin;

            warpedVol=vol;
            warpedKinVol=kinvol;
            warpedVol(:,:,ZNull)=[];
            warpedVol(YNull,:,:)=[];
            warpedVol(:,XNull,:)=[];
            warpedKinVol(:,:,ZNull)=[];
            warpedKinVol(YNull,:,:)=[];
            warpedKinVol(:,XNull,:)=[];

        otherwise
            error('unknown trans type');
    end

    %         case 'full'
    %
    %             warpedVol=imwarp(vol,inputRef,tform);
    %             warpedMaskedVol=imwarp(maskedVol,inputRef,tform);
    %
    %             warpedKinVol=imwarp(kinvol,inputRef,tform);
    %             warpedMaskedKinVol=imwarp(maskedKin,inputRef,tform);
    %     end


    %     imwrite(uint8(255*mat2gray(max(vol,[],3))),[outputDirDemo filesep  'vol_' num2str(fIdx,'%04d') '.png']) % Names the file stitched001..., in /stitched/
    %     imwrite(uint8(255*mat2gray(max(warpedVol,[],3))),[outputDirDemo filesep  'warpedVol_' num2str(fIdx,'%04d') '.png']) % Names the file stitched001..., in /stitched/
    %     imwrite(uint8(255*mat2gray(max(warpedMaskedVol,[],3))),[outputDirDemo filesep  'warpedMaskedVol_' num2str(fIdx,'%04d') '.png']) % Names the file stitched001..., in /stitched/
    %     imwrite(uint8(255*mat2gray(max(warpedKinVol,[],3))),[outputDirDemo filesep  'warpedKinVol_' num2str(fIdx,'%04d') '.png']) % Names the file stitched001..., in /stitched/
    %     imwrite(uint8(255*mat2gray(max(warpedMaskedKinVol,[],3))),[outputDirDemo filesep  'warpedKinMaskedVol_' num2str(fIdx,'%04d') '.png']) % Names the file stitched001..., in /stitched/


    %% Create MIPS for each channel, fuse mask and full volume
    switch p.transType
        case 'none'
        case 'transCrop'
            ZRatio=MD.pixelSizeZ_/MD.pixelSize_;
        otherwise
            ZRatio=1;
    end;

    % Create MIP of ROI and context
    [fullmaxXY,fullmaxZY,fullmaxZX,~]=computeMIPs(warpedVol,ZRatio, ...
        prctile(warpedVol(:),p.intMinPrctil(1)),prctile(warpedVol(:),p.intMaxPrctil(1)));
    [fullmaxXYKin,fullmaxZYKin,fullmaxZXKin,~]=computeMIPs(warpedKinVol,ZRatio, ...
        prctile(warpedKinVol(:),p.intMinPrctil(2)),prctile(warpedKinVol(:),p.intMaxPrctil(2)));

    [maxXY,maxZY,maxZX,~]=computeMIPs(warpedMaskedVol,ZRatio, ...
        prctile(warpedMaskedVol(:),p.intMinPrctil(1)),prctile(warpedMaskedVol(:),p.intMaxPrctil(1)));
    [maxXYKin,maxZYKin,maxZXKin,~]=computeMIPs(warpedMaskedKinVol,ZRatio, ...
        prctile(warpedMaskedKinVol(:),p.intMinPrctil(2)),prctile(warpedMaskedKinVol(:),p.intMaxPrctil(2)));

    % Fuse ROI and context
    maxXY(maxXY==0)=fullmaxXY(maxXY==0);
    maxZY(maxZY==0)=fullmaxZY(maxZY==0);
    maxZX(maxZX==0)=fullmaxZX(maxZX==0);
    maxXYKin(maxXYKin==0)=fullmaxXYKin(maxXYKin==0);
    maxZYKin(maxZYKin==0)=fullmaxZYKin(maxZYKin==0);
    maxZXKin(maxZXKin==0)=fullmaxZXKin(maxZXKin==0);

    %% Resize and fuse channel MIPS
    maxMIPSize=400;
    [sX,sY,sZ]=size(warpedMaskedVol);
    if(strcmp(p.transType,'none'))
        sZ=sZ*MD.pixelSizeZ_/MD.pixelSize_;
    end
    resizeScale=maxMIPSize/max([sX,sY,sZ]);

    rMax=imresize(maxXY,resizeScale);
    rmaxKin=imresize(maxXYKin,resizeScale);
%    mipSize=size(rMax);
    XYProj=renderChannel(rMax,rmaxKin,p.channelRender);

    rMax=imresize(maxZY,resizeScale);
    rmaxKin=imresize(maxZYKin,resizeScale);
    ZYProj=renderChannel(rMax,rmaxKin,p.channelRender);

    rMax=imresize(maxZX,resizeScale);
    rmaxKin=imresize(maxZXKin,resizeScale);
    ZXProj=renderChannel(rMax,rmaxKin,p.channelRender);

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
    
    

