function ROI = LeverFrameSegment_texture(I)

MIN_CELL_RADIUS_um = 2; % um
MAX_CELL_RADIUS_um = 50;

    chan=channels(idChannel);
    im = MicroscopeData.Reader('imageData',CONSTANTS.imageData, 'chanList',chan, ...
        'timeRange',[t t], 'outType','single','prompt',false);
    im=Segment.polyfix(im);
    
    resolution_um=CONSTANTS.imageData.PixelPhysicalSize(1); % um per pixel
    min_radius_pixels = MIN_CELL_RADIUS_um / resolution_um;
    min_area_pixels = min_radius_pixels^2 * pi;
    
    max_radius_pixels=MAX_CELL_RADIUS_um/resolution_um;
    max_area_pixels=max_radius_pixels^2*pi;
    
    se=strel('disk',2*ceil(min_radius_pixels));
    imf=entropyfilt(im,getnhood(se));
    % imf=ImProc.EntropyFilter(im,ones(51,51,1));
    imf=mat2gray(imf);
    
    % morphological gradient
    seBig=strel('disk',8*ceil(min_radius_pixels));
    seSmall=strel('disk',floor(min_radius_pixels));
    ig = imdilate(imf,seBig)-imerode(imf,seSmall);
    ig=mat2gray(ig);
    bwig=im2bw(ig,graythresh(ig));
    %
    bw=bwig;
    % imagesc(im);colormap(gray);hold on
    bwHO=imfill(bw,'holes') & ~bw;
    [Cells, bwMask]=extractCells(bwHO,'r', t, min_area_pixels/2, max_area_pixels,[],false);
    
    bwFill=bw;
    for i=1:2*min_radius_pixels
        %      imagesc(bwFill);drawnow;pause(1);title(num2str(i))
        bwFill = imdilate(bwFill,strel('disk',1));
        bwholes=(imfill(bwFill,'holes') & ~bwFill);
        bwHO = bwHO | bwholes;
        bwFill=bwFill|bwholes;
        [newCells bwMask]=extractCells(bwHO,'g', t, min_area_pixels, max_area_pixels,bwMask,false);
        Cells=[Cells newCells];
    end
    
    tElapsed=toc;
    fprintf(1,'segmented frame %d channel %d, found %d hulls, time=%f\n',t,chan,length(Cells),tElapsed);


%%
function [Cells, bwMask]=extractCells(bw,color, t, min_area_pixels, max_area_pixels,...
    bwMask,bDilate)

if nargin<6 || isempty(bwMask)
    bwMask=0*bw;
end

if nargin<7
    bDilate=false;
end

Cells=[];
% stats=regionprops(bw,'eccentricity','centroid');
[L, num]=bwlabel(bw);
for n=1:num
    
    [r c]=find(L==n);
    idx=sub2ind(size(bw),r,c);
    if length(r)>max_area_pixels
        bw(idx)=0;
        continue;
    end
    
    if length(r)<min_area_pixels
        bw(idx)=0;
        continue
    end
    
    ch = convhull(c,r);
    
    bwp=roipoly(bw,c(ch),r(ch));
    idxPoly=find(bwp);
    if any(bwMask(idxPoly))
        bw(idx)=0;
        continue
    end
    bwMask(idxPoly)=1;
    
    % do we want to dilate it?
    if bDilate
        bwn=0*bw;
        bwn(idx)=1;
        bwn=imdilate(bwn,strel('disk',15));
        [r, c]=find(bwn);
        ch = convhull(c,r);
    end
    
    %     plot(c,r,'.m');
    %       plot(c(ch),r(ch),'-','color',color,'linewidth',2);
    %      text(stats(n).Centroid(1),stats(n).Centroid(2),num2str(stats(n).Eccentricity),'color','w');
    newCell=[];
    newCell.time=t;
    newCell.centroid=[mean([c,r]) 0];
    newCell.surface=[c(ch),r(ch)];
    newCell.pts=[c,r];
    Cells=[Cells newCell];
    %     drawnow
end

