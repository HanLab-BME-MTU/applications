%% ---- Parameters --- %%

outDir = 'A:\Papers\3d bio paper\Figure Panels';


figExpProps = {'DPI',600};

saveFigs = true;%Save figs to disk. Disable so faster for development.



axProps = {'FontSize',15};

labProps = {'FontSize',15};

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% Segmentation, Skeletonization / Branch Det and Curvature Figures ---------- %%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

movPaths = {'L:\nih\data_2_2011_fast_and_c3\Faster\Faster3\Faster3a3\movieData.mat',...
            'L:\nih\TdTmCAAX7-bleb\Blebbistatin\18a\s4\movieData.mat'};
                
        
movNames = {'Control','Blebbistatin'};

views = [-106.5 30;
          95.5 40];

curvToMake = {'LAcurv','LAgauss','LAmean'};
curvNames = {'Maximum Absolute PC','Gaussian','Mean'};
nCurv = numel(curvToMake);

nMov = numel(movPaths);
%%

for iMov = 1:nMov

    %% --- Load mask for these figures -- - %

    iChan =1;
    iFrame = 1;

    MD = MovieData.load(movPaths{iMov},0);

    pixXY = MD.pixelSize_;
    pixZ = MD.zSpacing_;


    iSegProc = MD.getProcessIndex('SegmentationProcess3D',1,1);
    maskFiles = MD.processes_{iSegProc}.getOutMaskFileNames(iChan);
    mask{iMov} = tif3Dread([MD.processes_{iSegProc}.outFilePaths_{iChan} filesep maskFiles{1}{iFrame}]);      

    mask{iMov} = make3DImageVoxelsSymmetric(mask{iMov},pixXY,pixZ);


    %% --- Get curvature data for these figures -- - %

    %Just-re run this here so we get the locally-averaged curvature
    %values.(these cells were analyzed before these were saved to disk)
    sampRadPix = 2e3 / pixXY;

    tic,maskProp{iMov} = analyze3DMaskGeometry(mask{iMov},.15,[],sampRadPix);toc
    %Swap X and Y so it matches with display in imaris
    maskProp{iMov}.SmoothedSurface.vertices = maskProp{iMov}.SmoothedSurface.vertices(:,[2 1 3]);

    %Convert relevant curvature measures to microns
    maskProp{iMov}.locAvgCurv.LocMeanMaxAbsCurvature = maskProp{iMov}.locAvgCurv.LocMeanMaxAbsCurvature ./ (pixXY/1e3);
    maskProp{iMov}.locAvgCurv.LocMeanGaussianCurvature = maskProp{iMov}.locAvgCurv.LocMeanGaussianCurvature ./ (pixXY/1e3)^2;
    maskProp{iMov}.locAvgCurv.LocMeanMeanCurvature = maskProp{iMov}.locAvgCurv.LocMeanMeanCurvature ./ -(pixXY/1e3);
    



    %% --------- Curvature overlay figures ------- %%

    
    for j = 1:nCurv

        currFig = fsFigure(.75);

        showMaskSurfaceProp(maskProp{iMov},curvToMake{j});

        axHan = get(currFig,'CurrentAxes');
        
        if iMov == 1
            cRanges(j,:) = caxis;
        else
            caxis(cRanges(j,:))
        end
        

        colormap(jet(1024))
        colorbar
        allChild = get(axHan,'Children');
        lh1 = allChild(strcmp('light',get(allChild,'type')));
        lh2 = light;
        set(lh2,'Position',[-1 1 1])
        lighting phong


        view(views(iMov,1),views(iMov,2))

        axis off
        set(currFig,'Color','w')


        set(axHan,axProps{:})

        figName = [outDir filesep 'curvature colormapped surface ' movNames{iMov} ' ' curvNames{j}];
        if saveFigs,mfFigureExport(currFig,figName,figExpProps{:});end

    end
    
end








%% --------- Example Geometries for Curvature Illustration -------- %%



w = 100;
delta = 1;
r = sqrt(2*w^2)*1.1;%So we don't reach the singularities in our stupid coord system
x = -w:delta:w;
[X,Y] = meshgrid(x);



curvFun = {@(x,y)(sqrt(r^2 - x .^2 - y .^2)),
           @(x,y)(-sqrt(r^2 - x .^2 - y .^2)),           
           @(x,y)(sqrt(r^2 - x .^2))};           
curvCatNames = {'k1 pos k2 pos',
            'k2 neg k1 neg',
            'k1 pos k2 0',
            'k1 pos k2 neg'};
curvCatColors = [0 0 1 ;
                 1 0 1 ;
                 0 1 0 ;
                 1 1 0 ];                 
        
       
nCurvCat = numel(curvFun);

for iFun = 1:3%Because for the saddle we use polar coord
    
    currFig = fsFigure(.5);
    
    z = curvFun{iFun}(X,Y);
    
    %z(R > 50) = NaN;
    %Z(R > 50) = NaN;
    
    sHan = surf(X,Y,z,'EdgeColor','none','FaceColor',curvCatColors(iFun,:),'FaceAlpha',.7);
    hold on
    plot3(x,zeros(2*w+1,1),z(w+1,:),'k','LineWidth',3)
    plot3(zeros(2*w+1,1),x,z(:,w+1),'w','LineWidth',3)
    
    lighting phong
    
    lh = light;
    lh2 = light;
    set(lh2,'Position',[-1 1 1])
    
                    
    axis off
    axis equal
    set(currFig,'color','w')       
    
    figName = [outDir filesep 'surface for curv illustration ' curvCatNames{iFun} ];
    
    if saveFigs,mfFigureExport(currFig,figName,figExpProps{:});end
    
    
end


%And because I'm anal - the hyperbolic saddle above doesn't have constant
%curvature along the plotted principle curves, so we make this saddle,
%which does
iFun = 4;

R = sqrt(X .^2 + Y .^2);
T = atan(Y ./ X);%theta

locR = sec(2*T)*r;%Radius varius with angle, but is constant at any given angle
z = sign(locR) .* sqrt(locR .^ 2 - R .^2) - locR;%Signed radius.
z(w+1,w+1) = 0;%Fix the stupid discontinuity at origin

currFig = fsFigure(.5);
sHan = surf(X,Y,z,'EdgeColor','none','FaceColor',curvCatColors(4,:),'FaceAlpha',.6);
hold on
plot3(x,zeros(2*w+1,1),z(w+1,:),'k','LineWidth',3)
plot3(zeros(2*w+1,1),x,z(:,w+1),'w','LineWidth',3)

lighting phong

lh = light;
lh2 = light;
set(lh2,'Position',[-1 1 1])

axis off
axis equal
set(currFig,'color','w')       
    

figName = [outDir filesep 'surface for curv illustration ' curvCatNames{iFun} ];
    
if saveFigs,mfFigureExport(currFig,figName,figExpProps{:});end
    


%% ----- Curvature Category overlay + closeups ------- %%

% ---------- Setup the categories and colors ----- %



for iMov = 1:nMov

    
    K = maskProp{iMov}.locAvgCurv.LocMeanGaussianCurvature;
    nFace = numel(K);   
    H = maskProp{iMov}.locAvgCurv.LocMeanMeanCurvature;
    
    if iMov == 1;
        %So we use the same limits for all movies
        hSmall = std(H);
        kSmall = std(K);
    end

    curvCatKLims = [ kSmall Inf;
                    -kSmall Inf;
                    -kSmall kSmall;
                    -Inf -kSmall];
    curvCatHLims = [hSmall  Inf;
                    -Inf    0;
                     0      Inf;
                    -hSmall Inf];

    curvCatFuns = {@(x)(K > kSmall & H > hSmall),...
                   @(x)(K > -kSmall & H < 0),...
                   @(x)(K > -kSmall & K < kSmall & H > 0),...
                   @(x)(K < -kSmall & H > -hSmall & H < Inf)};           

    nCurvCat = numel(curvCatFuns);
    ptCurvCat = nan(nFace,nCurvCat);
    perPointCols = zeros(nFace,3);

    
    
    currFig = figure;
    hold on
    for iCat = 1:nCurvCat    
        currPts = K > curvCatKLims(iCat,1) & K < curvCatKLims(iCat,2) & H > curvCatHLims(iCat,1) & H < curvCatHLims(iCat,2);

        ptCurvCat(currPts,iCat) = iCat;%Do it this way so we can check for category overlap

        plot(H(currPts),K(currPts),'.','color',curvCatColors(iCat,:));

        perPointCols(currPts,:) = repmat(curvCatColors(iCat,:),[nnz(currPts) 1]);
    end
    ptNoCat = sum(~isnan(ptCurvCat),2) == 0;

    %Do local-averaging of colors
    perPointColsAvg = nan(nFace,3);
    nVert = size(maskProp{iMov}.SmoothedSurface.vertices,1);
    faceCent = squeeze(mean(reshape(maskProp{iMov}.SmoothedSurface.vertices(maskProp{iMov}.SmoothedSurface.faces(:),:),[nFace,3,3]),2));%Get barycenter of each face
    tic
    allNeighb = KDTreeBallQuery(faceCent,faceCent,2.5);%Get faces close to each face
    toc
    tic
    for j = 1:nFace

        %Find all faces that share these vertices    
        perPointColsAvg(j,:) = mean(perPointCols(allNeighb{j},:),1);

    end
    toc




    catStr = arrayfun(@(x)([ 'Category ' num2str(x)]),1:nCurvCat,'Unif',0);
    catStr = [ catStr {'No Category'}];
    plot(H(ptNoCat),K(ptNoCat),'.k')
    plot([0 0],ylim,'--k')
    plot(xlim,[0 0],'--k')
    legend(catStr)
    xlabel('Mean Curvature')
    ylabel('Gaussian Curvature')
    figName=  [outDir filesep 'curvature categories scatter plot ' movNames{iMov}];
    if saveFigs,mfFigureExport(currFig,figName,figExpProps{:});end


    disp([num2str(nFace)  ' points total'])
    disp([num2str(nnz(sum(double(~isnan(ptCurvCat)),2)>1)) ' points in more than one category'])
    disp([num2str(nnz(sum(double(~isnan(ptCurvCat)),2)==0)) ' points in zero categories'])

    currFig = figure;
    tmp =ptCurvCat;
    tmp(isnan(tmp)) = 0;
    tmp = max(tmp,[],2);
    hist(tmp(:),0:nCurvCat)
    xlabel('Category # (0 is uncategorized)')
    ylabel('# Surface vertices')
    figName=  [outDir filesep 'curvature categories histogram ' movNames{iMov}];
    if saveFigs,mfFigureExport(currFig,figName,figExpProps{:});end
    
    %% ----- Make 2D histogram figures ---- %%
    
    currFig = figure;
    nBins = 100;
    if iMov == 1
      xl = [-.7 1.37];
        yl = [-1 1.5];
        C{1} = linspace(xl(1),xl(2),nBins);
        C{2} = linspace(yl(1),yl(2),nBins);
        
    end
    N = hist3([H K],C);
    
    imHan = imagesc(C{1},C{2},log10(N'+1));
    hold on
    axHan = get(currFig,'CurrentAxes');
    if iMov == 1
        %xl = prctile(H,[.05 99.95]);
        %xl = [min(H(:)) max(H(:))];
        %yl = prctile(K,[.05 99.95]);
        %yl = [min(K(:)) max(K(:))];
      
    end
    xlim(xl);
    ylim(yl);
    axis xy
    plot(xlim,[0 0],'--r')
    plot([0 0],ylim,'--r')
    saturateImageColormap(axHan,1);
    colormap gray
    xlabel('Mean Curvature, 1/microns')
    ylabel('Gaussian Curvature, 1/microns^2')
    colorbar
    
    figName=  [outDir filesep 'mean vs gaussian curvature 2D histogram ' movNames{iMov}];
    if saveFigs,mfFigureExport(currFig,figName,figExpProps{:});end
    
    %Now overlay the categories
    for iCat = 1:nCurvCat
        
        currX = [curvCatHLims(iCat,1) curvCatHLims(iCat,2) curvCatHLims(iCat,2) curvCatHLims(iCat,1)];
        currX(currX > xl(2)) = xl(2);
        currX(currX < xl(1)) = xl(1);
        currY = [curvCatKLims(iCat,1) curvCatKLims(iCat,1) curvCatKLims(iCat,2) curvCatKLims(iCat,2)];
        currY(currY > yl(2)) = yl(2);
        currY(currY < yl(1)) = yl(1);

        fill(currX,currY,curvCatColors(iCat,:),'FaceAlpha',.15,'EdgeColor',curvCatColors(iCat,:))
    
    end

     figName=  [outDir filesep 'mean vs gaussian curvature 2D histogram ' movNames{iMov} ' with categories'];
    if saveFigs,mfFigureExport(currFig,figName,figExpProps{:});end
    
    %% ----- Make the geometry overlay figures ----- %
    
    currFig = fsFigure(.75);
    patchHan = patch(maskProp{iMov}.SmoothedSurface,'FaceColor','flat','EdgeColor','none','FaceVertexCData',perPointColsAvg,'VertexNormals',maskProp{iMov}.SurfaceNorms(:,[2 1 3]));
    axis equal

    axHan = get(currFig,'CurrentAxes');
    allChild = get(axHan,'Children');
    lh1 = allChild(strcmp('light',get(allChild,'type')));
    lh2 = light;
    set(lh2,'Position',[-1 1 1])
    lighting phong


    view(views(iMov,1),views(iMov,2))

    axis off
    set(currFig,'Color','w')

    set(axHan,axProps{:})


    figName=  [outDir filesep 'curvature categories surface overlay whole cell ' movNames{iMov}];
    if saveFigs,mfFigureExport(currFig,figName,figExpProps{:});end

    %% ---- And the closeups ----- %

    closeCamPos = 1e3*[-1.3188 2.3699 1.138;
                      -1.3741    2.4513    1.1232];
    closeCamAng = [.9406;
                    1.7915];

    closeCamTarg = [ 271.8201  140.7660   33.6317;
                    125.0513  291.1307   20.2313];


    nClose = size(closeCamPos,1);

    for iClose = 1:nClose

        set(axHan,'CameraPosition',closeCamPos(iClose,:))
        set(axHan,'CameraViewAngle',closeCamAng(iClose))
        set(axHan,'CameraTarget',closeCamTarg(iClose,:));

        figName=  [outDir filesep 'curvature categories surface overlay ' movNames{iMov} ' closeup ' num2str(iClose)];
        if saveFigs,mfFigureExport(currFig,figName,figExpProps{:});end
    end
    
end


%% ------------- Branch Varation around a mean bullshit figures ------ %%

skPPFile = load('L:\nih\Post Processing Myo Inhib and WT\Skeleton Post Proc\WT\pruned skeleton post processing.mat');

currFig = fsFigure(.5);
hold on

nMovies = numel(skPPFile.nTipsMeanPerMovThresh);

%iMovShow = [ 1 3 17];%4
iMovShow = [ 3 4 8];%13
%iMovShow = 1:nMovies

nShow = numel(iMovShow);
movColors = lines(nMovies);

allMaxT = cellfun(@max,skPPFile.allTData);

maxT = max(allMaxT(iMovShow));

minMaxT = min(allMaxT(iMovShow));


for iMov = iMovShow
    %Do the mean/std bands per movie first so they're in the background
    mP = skPPFile.nTipsMeanPerMovThresh(iMov) + skPPFile.nTipsSTDPerMovThresh(iMov);
    mM = skPPFile.nTipsMeanPerMovThresh(iMov) - skPPFile.nTipsSTDPerMovThresh(iMov);
    fill([0 maxT maxT 0],[mP mP mM mM],movColors(iMov,:),'FaceAlpha',.2,'EdgeColor','none');
    
end

legend(arrayfun(@(x)(['Cell ' num2str(x)]),1:nShow,'Unif',0))

for iMov = iMovShow
     
     
    plot(skPPFile.allTData{iMov},skPPFile.nTipsPerFrameThresh{iMov},'-','Color',movColors(iMov,:),'LineWidth',3)
            
    %nTipsMeanPerMovThresh(iMov) = mean(nTipsPerFrameThresh{iMov});
    %nTipsSTDPerMovThresh(iMov) = std(nTipsPerFrameThresh{iMov});            
    
end

xlabel('Time, Seconds',labProps{:})
ylabel('Branch Number',labProps{:})    
axHan = get(currFig,'CurrentAxes');
set(axHan,axProps{:})

xlim([0 minMaxT])

figName=  [outDir filesep 'example branch count over time and STD'];
if saveFigs,mfFigureExport(currFig,figName,figExpProps{:});end

%% ----- Branch variation histogram figures ----- %%

currFig = figure;

%Get per-movie % coefficient of variation

[n,c] = hist(skPPFile.cvPctPerMovnTipsThresh,10);
bar(c,n,'b','EdgeColor','b')
hold on
xlabel('% Coefficient of Variation',labProps{:})
ylabel('# Cells',labProps{:})
title({'CV%, Thresholded Branch Tip #, per cell average/STD over all frames',...
    ['n = ' num2str(nMovies) ' cells'],...
    ['Mean = ' num2str(mean(skPPFile.cvPctPerMovnTipsThresh))]});
%^    ['Population CV% of per-cell means = ' num2str(cvPctAllMovnTipsThresh) '%']});
ylim([0 max(ylim)*1.1])
xlim([0 100])

axHan = get(currFig,'CurrentAxes');
set(axHan,axProps{:})


figName = [outDir filesep 'branch count percent CV per cell over time '];
if saveFigs,mfFigureExport(currFig,figName,figExpProps{:});end


currFig = figure;
hold on
[n,c] = hist(skPPFile.nTipsMeanPerMovThresh,10);
bar(c,n,'b','EdgeColor','b')
xlabel('Mean Branch Number, All Timepoints',labProps{:})
ylabel('# Cells',labProps{:})
title(['CV% Mean Branch Tip # Across cells = ' num2str(cvPctAllMovnTipsThresh) '%'])
ylim([0 max(ylim)*1.1])

axHan = get(currFig,'CurrentAxes');
set(axHan,axProps{:})

figName = [outDir filesep 'branch count mean per cell all frames and CV across cells'];
if saveFigs,mfFigureExport(currFig,figName,figExpProps{:});end


%% ----------- Whole-Condition curvature averages bar-graphs

parDir = 'L:\nih\Post Processing Myo Inhib and WT\Mask Geometry';
fName = 'combined stats and hists all movies.mat';


condNames = {'WT','Blebbistatin'};
nCond = numel(condNames);
condCols = [0 0 1; 1 0 0 ];

curvVars = {'gcMeanPerMov','mcMedPerMov','mapcMedPerMov'};

curvNames = {'Median Gaussian Curvature','Median Mean Curvature','Median Maximum Absolute Curvature'};
nCurv = numel(curvNames);

%We convert to microns so the nubmers are more reasonable, and flip the
%sign on the mean curvature
curvConv = [1e3^2 -1e3 1e3];
curvUnits = {'microns^-^2','microns^-^1','microns^-^1'};

combCurvData = cell(nCond,1);
for j = 1:nCond
    combCurvData{j} = load([parDir filesep condNames{j} filesep fName]);
end
   
%%



for iCurv = 1:nCurv        
    
    clear currMeanPerCell
    
    currFig = figure;
    hold on

    nCellStr = {};
    
    for iCon = 1:nCond
        
        currMeanPerCell{iCon} = cellfun(@mean,combCurvData{iCon}.(curvVars{iCurv})) .* curvConv(iCurv);
        
        if iCon == 1
            %First two cells are from oldest data and have wrong pixel size, we just exclude them
            %for now.
            currMeanPerCell{iCon} = currMeanPerCell{iCon}(3:end);
        end
                
        allMean(iCurv,iCon) = mean(currMeanPerCell{iCon});
        nCellsPer(iCon) = numel(currMeanPerCell{iCon});
        
        nCellStr{iCon} =  [condNames{iCon} ' n=' num2str( nCellsPer(iCon))];
        
        semPerCurvCond(iCurv,iCon) = std(currMeanPerCell{iCon}) / sqrt(nCellsPer(iCon));
    
        barHan(iCurv,iCon) = bar(iCon,allMean(iCurv,iCon),.9,'FaceColor',condCols(iCon,:));
                
    end
    
    yl = [min(allMean(iCurv,:)) max(allMean(iCurv,:))];
    if iCurv ~= 1
        yl = yl .* [.95 1.05] ;
    else
        yl = yl .* [1.5 .7]; %Because the gaussian combined stats are negative, and have large SEM
    end
    
    [h,p] = ttest2(currMeanPerCell{1},currMeanPerCell{2},.05,'both','unequal');
    
    if h
        sigStr = '*';
    else
        sigStr = 'N.S.';
    end
    
    ylim(yl)
    
    title({['Mean of Per-Cell ' curvNames{iCurv}],...
        'Bars show +/- SEM',...
        [sigStr ', p=' num2str(p)] nCellStr{:}});
    %legend(condNames,'Location','NorthEastOutside')        
    ylabel(curvUnits{iCurv})
    axHan = get(currFig,'CurrentAxes');
    set(axHan,'XTick',1:nCond)
    set(axHan,'XTickLabel',condNames)
    
    %DO last so legend works
    for iCon = 1:nCond
        eBarHan(iCurv,iCon) = errorbar(iCon,allMean(iCurv,iCon),semPerCurvCond(iCurv,iCon),'k');
    end
    
    figName = [outDir filesep 'bar graph comparison whole cell combined ' curvNames{iCurv}];
    if saveFigs,mfFigureExport(currFig,figName,figExpProps{:});end
    
end










