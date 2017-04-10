function projImages=project1D(MD,dynPoligonISO,varargin)
% tracks are in original FoF and projected  in the manifold
% described by polygon.
% Poligon described the projected space (without dilated boundaries), It is composed of Tracks.
% in 1D, it is only composed of two tracks.

ip = inputParser;
ip.CaseSensitive = false;
ip.KeepUnmatched = true;
ip.addOptional('dynPoligonREF',dynPoligonISO);
ip.addOptional('tracks',[]);
ip.addOptional('crop','manifold');
ip.addOptional('transType','affineOnePass');
ip.addOptional('FoF',[]);
ip.addOptional('channelRender','greenRed');
ip.addOptional('name',[]);
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
inputRef=imref3d([ MD.getDimensions('Y') MD.getDimensions('X') MD.getDimensions('Z')], ...
    [1 MD.getDimensions('X')],[1 MD.getDimensions('Y')],[1 MD.getDimensions('Z')*MD.pixelSizeZ_/MD.pixelSize_]);


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
        OrigSubFrame=orig(dynPoligonREF(iP).f,:);
        minX=floor(min(min(dynPoligonREF(iP).x),minX));
        minY=floor(min(min(dynPoligonREF(iP).y),minY));
        minZ=floor(min(min(dynPoligonREF(iP).z),minZ));
        maxX=ceil(max(max(dynPoligonREF(iP).x),maxX));
        maxY=ceil(max(max(dynPoligonREF(iP).y),maxY));
        maxZ=ceil(max(max(dynPoligonREF(iP).z),maxZ));
    end
    % minZ=floor(minZ*MD.pixelSize_/MD.pixelSizeZ_);
    % maxZ=ceil(maxZ*MD.pixelSize_/MD.pixelSizeZ_);

    maxXBorder=min(maxX+cubeHalfWidth,MD.getDimensions('X'));
    maxYBorder=min(maxY+cubeHalfWidth,MD.getDimensions('Y'));
    maxZBorder=min(maxZ+cubeHalfWidth,(MD.pixelSizeZ_/MD.pixelSize_)*MD.getDimensions('Z'));
    minXBorder=max(minX-cubeHalfWidth,1);
    minYBorder=max(minY-cubeHalfWidth,1);
    minZBorder=max(minZ-cubeHalfWidth,1);

    maxXBorder=(maxX+cubeHalfWidth);
    maxYBorder=(maxY+cubeHalfWidth);
    maxZBorder=(maxZ+cubeHalfWidth);
    minXBorder=(minX-cubeHalfWidth);
    minYBorder=(minY-cubeHalfWidth);
    minZBorder=(minZ-cubeHalfWidth);
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
outputDirDemo=[MD.outputDirectory_ filesep '1DProjection' filesep p.name filesep 'volDemo' ];

system(['mkdir  -p ' outputDirSlices1]);
system(['mkdir  -p ' outputDirDemo]);
for fIdx=1:MD.nFrames_
    vol=MD.getChannel(1).loadStack(fIdx);
    kinvol=MD.getChannel(2).loadStack(fIdx);
    % Collect relative frameIdx
    pIndices=nan(1,length(dynPoligonISO));
    for polIdx=1:length(dynPoligonISO)
        F=dynPoligonISO(polIdx).f;
        pIdx=find(F==fIdx);
        if isempty(pIdx)
            if(fIdx>max(F))   pIdx=max(F);  else   pIdx=min(F); end;
        end
        pIndices(polIdx)=pIdx;
    end;
    PCurrent=[dynPoligonISO(1).x(pIndices(1)) dynPoligonISO(1).y(pIndices(1)) dynPoligonISO(1).z(pIndices(1))];
    KCurrent=[dynPoligonISO(2).x(pIndices(2)) dynPoligonISO(2).y(pIndices(2)) dynPoligonISO(2).z(pIndices(2))];

    % Building mask for both channel on the whole volume
    % NOTE: Could masking use imwarp for speed ?
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
%     outputDir=[outputDirDemo filesep  'mask'];mkdir(outputDir);
%     imwrite(uint8(255*mat2gray(max(maskedVol,[],3))),[outputDir filesep  'mask_' num2str(fIdx,'%04d') '.png']) % Names the file stitched001..., in /stitched/

    maskedKin=kinvol;
    maskedKin(~mask)=0;
%     outputDir=[outputDirDemo filesep  'maskKin_'];mkdir(outputDir);
%     imwrite(uint8(255*mat2gray(max(maskedKin,[],3))),[outputDir filesep  'maskKin_' num2str(fIdx,'%04d') '.png']) % Names the file stitched001..., in /stitched/


    % If needed the map must rotated before cropped (efficiency)
    % Rotation depends on the FrameOfRef associated to the tracks the compose the dynanimc polygon
    % Cropping area according to the polygon OVER TIME plus added vizualiation margin
    % Rotation will use imwarp
    % Can we use imwar for cropping too ?

    tform=affine3d();

    if(~isempty(p.FoF))
        B=p.FoF.getBase(fIdx);
        tform.T(4,[1 2 3])=(-p.FoF.origin(fIdx,:)+p.FoF.origin(1,:))*B;
        tform.T(1:3,1:3)=B;
        %
        tformTransOnly=affine3d();
        tformTransOnly.T(4,[1 2 3])=(-p.FoF.origin(fIdx,:));

                %
        tformRelTransOnly=affine3d();
        tformRelTransOnly.T(4,[1 2 3])=(-p.FoF.origin(1,:)+p.FoF.origin(fIdx,:));

        tformRotOnly=affine3d();
        B=p.FoF.getBase(fIdx);
        tformRotOnly.T(1:3,1:3)=B;

        tformRotOnlyInit=affine3d();
        B=p.FoF.getBase(1);
        tformRotOnlyInit.T(1:3,1:3)=B;
    end

    switch p.transType
        case 'affine'
            maxXBorderFull=MD.getDimensions('X');
            maxYBorderFull=MD.getDimensions('Y');
            maxZBorderFull=MD.getDimensions('Z')*(MD.pixelSizeZ_/MD.pixelSize_);
            minXBorderFull=1;
            minYBorderFull=1;
            minZBorderFull=1;

            border=300;
            minXBorderCurr=minXBorderFull - p.FoF.origin(fIdx,1)-border; maxXBorderCurr=maxXBorderFull - p.FoF.origin(fIdx,1)+border;
            minYBorderCurr=minYBorderFull - p.FoF.origin(fIdx,2)-border; maxYBorderCurr=maxYBorderFull - p.FoF.origin(fIdx,2)+border;
            minZBorderCurr=minZBorderFull - p.FoF.origin(fIdx,3)-border; maxZBorderCurr=maxZBorderFull - p.FoF.origin(fIdx,3)+border;

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
            minXBorderCurr=minXBorder ;%+ p.FoF.origin(fIdx,1) - p.FoF.origin(1,1);
            maxXBorderCurr=maxXBorder ;%+ p.FoF.origin(fIdx,1) - p.FoF.origin(1,1);
            minYBorderCurr=minYBorder ;%+ p.FoF.origin(fIdx,2) - p.FoF.origin(1,2);
            maxYBorderCurr=maxYBorder ;%+ p.FoF.origin(fIdx,2) - p.FoF.origin(1,2);
            minZBorderCurr=minZBorder ;%+ p.FoF.origin(fIdx,3) - p.FoF.origin(1,3);
            maxZBorderCurr=maxZBorder ;%+ p.FoF.origin(fIdx,3) - p.FoF.origin(1,3);

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

            minXBorderCurr=minXBorderFull - p.FoF.origin(fIdx,1); maxXBorderCurr=maxXBorderFull -  p.FoF.origin(fIdx,1);
            minYBorderCurr=minYBorderFull - p.FoF.origin(fIdx,2); maxYBorderCurr=maxYBorderFull  - p.FoF.origin(fIdx,2);
            minZBorderCurr=minZBorderFull - p.FoF.origin(fIdx,3); maxZBorderCurr=maxZBorderFull -  p.FoF.origin(fIdx,3);

            inputRef=imref3d([ MD.getDimensions('Y') MD.getDimensions('X') MD.getDimensions('Z')], ...
            [minXBorderCurr maxXBorderCurr], [minYBorderCurr maxYBorderCurr], [minZBorderCurr maxZBorderCurr]);

            disp(num2str(fIdx))
            minXBorderCurr=minXBorder ;%+ p.FoF.origin(fIdx,1) - p.FoF.origin(1,1);
            maxXBorderCurr=maxXBorder ;%+ p.FoF.origin(fIdx,1) - p.FoF.origin(1,1);
            minYBorderCurr=minYBorder ;%+ p.FoF.origin(fIdx,2) - p.FoF.origin(1,2);
            maxYBorderCurr=maxYBorder ;%+ p.FoF.origin(fIdx,2) - p.FoF.origin(1,2);
            minZBorderCurr=minZBorder ;%+ p.FoF.origin(fIdx,3) - p.FoF.origin(1,3);
            maxZBorderCurr=maxZBorder ;%+ p.FoF.origin(fIdx,3) - p.FoF.origin(1,3);

            % border=0;
            % minXBorderCurr=minXBorder + p.FoF.origin(fIdx,1) - p.FoF.origin(1,1)-border;
            % maxXBorderCurr=maxXBorder + p.FoF.origin(fIdx,1) -  p.FoF.origin(1,1)-border;
            % minYBorderCurr=minYBorder + p.FoF.origin(fIdx,2) - p.FoF.origin(1,2)-border;
            % maxYBorderCurr=maxYBorder + p.FoF.origin(fIdx,2) - p.FoF.origin(1,2)-border;
            % minZBorderCurr=minZBorder + p.FoF.origin(fIdx,3) - p.FoF.origin(1,3)-border;
            % maxZBorderCurr=maxZBorder + p.FoF.origin(fIdx,3) -  p.FoF.origin(1,3)-border;

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
            warpedMaskedVol=imwarp(maskedVol,inputRef,tformRotOnly,'OutputView',rotOutputRef);

            warpedKinVol=imwarp(kinvol,inputRef,tformRotOnly,'OutputView',rotOutputRef);
            warpedMaskedKinVol=imwarp(maskedKin,inputRef,tformRotOnly,'OutputView',rotOutputRef);
        case 'translation'
            disp(num2str(fIdx))
            minXBorderCurr=minXBorder ;%+ p.FoF.origin(fIdx,1) - p.FoF.origin(1,1);
            maxXBorderCurr=maxXBorder ;%+ p.FoF.origin(fIdx,1) - p.FoF.origin(1,1);
            minYBorderCurr=minYBorder ;%+ p.FoF.origin(fIdx,2) - p.FoF.origin(1,2);
            maxYBorderCurr=maxYBorder ;%+ p.FoF.origin(fIdx,2) - p.FoF.origin(1,2);
            minZBorderCurr=minZBorder ;%+ p.FoF.origin(fIdx,3) - p.FoF.origin(1,3);
            maxZBorderCurr=maxZBorder ;%+ p.FoF.origin(fIdx,3) - p.FoF.origin(1,3);

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

            minYBorderCurr=minYBorder + p.FoF.origin(fIdx,2); % canonical ref projection
            YNull= zeros(1,size(maskcrop,1));
            YNull(1:minYBorderCurr)=1;
            YNull(minYBorderCurr:end)=1;
            YNull=logical(YNull);

            minXBorderCurr=minXBorder + p.FoF.origin(fIdx,1); % canonical ref projection
            XNull= zeros(1,size(maskcrop,2));
            XNull(1:minXBorderCurr)=1;
            XNull(minXBorderCurr:end)=1;
            XNull=logical(XNull);

            minZBorderCurr=minZBorder + p.FoF.origin(fIdx,3); % canonical ref projection
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

    % normalize intensities


    [fullmaxXY,fullmaxZY,fullmaxZX,~]=computeMIPs(warpedVol,ZRatio,prctile(warpedVol(:),1),prctile(warpedVol(:),100));
    [fullmaxXYKin,fullmaxZYKin,fullmaxZXKin,~]=computeMIPs(warpedKinVol,ZRatio,prctile(warpedKinVol(:),99.9),prctile(warpedKinVol(:),100));

    [maxXY,maxZY,maxZX,~]=computeMIPs(warpedMaskedVol,ZRatio,prctile(warpedMaskedVol(:),1),prctile(warpedMaskedVol(:),100));
    [maxXYKin,maxZYKin,maxZXKin,~]=computeMIPs(warpedMaskedKinVol,ZRatio,prctile(warpedMaskedKinVol(:),99.9),prctile(warpedMaskedKinVol(:),100));

    maxXY(maxXY==0)=fullmaxXY(maxXY==0);
    maxZY(maxZY==0)=fullmaxZY(maxZY==0);
    maxZX(maxZX==0)=fullmaxZX(maxZX==0);
    maxXYKin(maxXYKin==0)=fullmaxXYKin(maxXYKin==0);
    maxZYKin(maxZYKin==0)=fullmaxZYKin(maxZYKin==0);
    maxZXKin(maxZXKin==0)=fullmaxZXKin(maxZXKin==0);


    %% Fuse channel MIPS add track and build mosaic.
    maxMIPSize=400;
    [sX,sY,sZ]=size(warpedMaskedVol);
    if(strcmp(p.transType,'none'))
        sZ=sZ*MD.pixelSizeZ_/MD.pixelSize_;
    end
    resizeScale=maxMIPSize/max([sX,sY,sZ]);


    rMax=imresize(maxXY,resizeScale);
    rmaxKin=imresize(maxXYKin,resizeScale);
    mipSize=size(rMax);
%     mipSize=[400 400];
%     rMaxXY=imresize(maxXY,mipSize);
%     rmaxXYKin=imresize(maxXYKin,mipSize);

    RGBThree=renderChannel(rMax,rmaxKin,p.channelRender);

    
    %% select tracks that starts or end in the manifold
    inMaskTrack=ones(1,length(tracks));
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
    %%
    tracksInMask=tracks(logical(inMaskTrack));
    
    %%
    myColormap=[[0 0 255]; [0 255 00]];
    tracksColors=uint8(myColormap);

    if(~isempty(tracksInMask))
      tracksXY=trackBinaryOverlay(RGBThree,[minXBorder maxXBorder],[minYBorder maxYBorder],tracksInMask,fIdx,ones(1,length(tracksInMask)),tracksColors);
    else
      tracksXY=RGBThree;
    end

    rMax=imresize(maxZY,resizeScale);
    rmaxKin=imresize(maxZYKin,resizeScale);

    RGBThree=renderChannel(rMax,rmaxKin,p.channelRender);

    if(~isempty(tracksInMask))
      capturedEB3ZY=tracksInMask.copy();
      for ebIdx=1:length(capturedEB3ZY)
        capturedEB3ZY(ebIdx).x=tracksInMask(ebIdx).z ;%*MD.pixelSize_/MD.pixelSizeZ_;
      end
      tracksColors=uint8(myColormap);
      tracksZY=trackBinaryOverlay(RGBThree,[minZBorder maxZBorder],[minYBorder maxYBorder],capturedEB3ZY,fIdx,ones(1,length(tracksInMask)),tracksColors);
    else
      tracksZY=RGBThree;
    end;

    tracksZY=permute(tracksZY,[2 1 3]);
    stripeSize=4;
    threeTop = [tracksXY, zeros(size(tracksXY,1), stripeSize,3), zeros(size(tracksXY,1), size(tracksZY,2),3)];


    rMax=imresize(maxZX,resizeScale);
    rmaxKin=imresize(maxZXKin,resizeScale);

    RGBThree=renderChannel(rMax,rmaxKin,p.channelRender);

    if(~isempty(tracks))
      capturedEB3ZX=capturedEB3ZY.copy();
      for ebIdx=1:length(capturedEB3ZX)
        capturedEB3ZX(ebIdx).y=tracksInMask(ebIdx).x;
      end
      tracksColors=uint8(myColormap);
      tracksZX=trackBinaryOverlay(RGBThree,[minZBorder maxZBorder],[minXBorder maxXBorder],capturedEB3ZX,fIdx,ones(1,length(tracksInMask)),tracksColors);
    else
      tracksZX=RGBThree;
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
