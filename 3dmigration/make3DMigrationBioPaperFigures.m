%% ---- Parameters --- %%

%outDir = 'A:\Papers\3d bio paper\Figure Panels';
outDir = 'W:\Hunter\orchestra_files_and_backup_merged\nih\Figures\Figure Panels New';


figExpProps = {'DPI',600};

saveFigs = true;%Save figs to disk. Disable so faster for development.



axProps = {'FontSize',15};

labProps = {'FontSize',15};

plotPars = {'LineWidth',2};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%% Segmentation, Skeletonization / Branch Det and Curvature Figures ---------- %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% movPaths = {'L:\nih\data_2_2011_fast_and_c3\Faster\Faster3\Faster3a3\movieData.mat',...
%             'L:\nih\TdTmCAAX7-bleb\Blebbistatin\18a\s4\movieData.mat'};

movPaths = {'W:\Hunter\orchestra_files_and_backup_merged\nih\data_2_2011_fast_and_c3\Faster\Faster3\Faster3a3\movieData.mat',...
            'W:\Hunter\orchestra_files_and_backup_merged\nih\TdTmCAAX7-bleb\Blebbistatin\18a\s4\movieData.mat'};
                
        
        
        
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

[curvCatKLims,curvCatHLims,curvCatNames,curvCatColors] = getCurvCategories;


w = 100;
delta = 1;
r = sqrt(2*w^2)*1.1;%So we don't reach the singularities in our stupid coord system
x = -w:delta:w;
[X,Y] = meshgrid(x);



curvFun = {@(x,y)(sqrt(r^2 - x .^2 - y .^2)),
           @(x,y)(-sqrt(r^2 - x .^2 - y .^2)),           
           @(x,y)(sqrt(r^2 - x .^2))};           
        
       
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
    
%     if iMov == 1;
%         %So we use the same limits for all movies
%         hSmall = std(H);
%         kSmall = std(K);
%     end
% 
%     curvCatKLims = [ kSmall Inf;
%                     -kSmall Inf;
%                     -kSmall kSmall;
%                     -Inf -kSmall];
%     curvCatHLims = [hSmall  Inf;
%                     -Inf    0;
%                      0      Inf;
%                     -hSmall Inf];
% 
%     curvCatFuns = {@(x)(K > kSmall & H > hSmall),...
%                    @(x)(K > -kSmall & H < 0),...
%                    @(x)(K > -kSmall & K < kSmall & H > 0),...
%                    @(x)(K < -kSmall & H > -hSmall & H < Inf)};           
   
               

    nCurvCat = numel(curvCatNames);
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


%% ---------- Curvature Category Histograms -------- %%

%Get limits and indices for curvature categories.
[curvCatKLims,curvCatHLims,curvCatNames,curvCatColors] = getCurvCategories;
nCurvCat = numel(curvCatNames);

ppFile = 'combined stats and hists all movies.mat';
ppDirs = {'W:\Hunter\orchestra_files_and_backup_merged\nih\Post Processing Myo Inhib and WT\Mask Geometry\WT',...
          'W:\Hunter\orchestra_files_and_backup_merged\nih\Post Processing Myo Inhib and WT\Mask Geometry\Blebbistatin';};
conNames = {'Untreated','Blebbistatin'};      
nCon = numel(conNames);
conColors = [ 0 0 0 ; 1 0 0 ];
conStyles = {'-','--'};
catField = 'meanFracPerMov';
nMovPer = nan(nCon,1);
isNorm = nan(nCon,nCurvCat);

testAlpha = .05;

for j = 1:nCon
    mgPP(j) = load([ppDirs{j} filesep ppFile]);
    nMovPer(j) = size(mgPP(j).(catField),1);    
    
    for k = 1:nCurvCat
        %Test for normality
        [isNorm(j,k),adPval(j,k)] = adtest(mgPP(j).(catField)(:,k));                
    end
    meanPer(j,:) = mean(mgPP(j).(catField),1);
    stdPer(j,:) = std(mgPP(j).(catField),[],1);
    semPer(j,:) = stdPer(j,:) ./ sqrt(nMovPer(j));
end

for k = 1:nCurvCat
    
    [pMW(k),hMW(k)] = ranksum(mgPP(1).(catField)(:,k),mgPP(2).(catField)(:,k),'Alpha',testAlpha);

end
isNorm = ~isNorm;

currFig = fsFigure(.5);
hold on
currAx = get(currFig,'CurrentAxes');

xEq = @(j,k)(j-1/6+(k-1)/3);

for j = 1:nCurvCat
    
    for k = 1:nCon
        currX = xEq(j,k);
        bHan = bar(currX,meanPer(k,j),.30);
        
        if j == 1 && k == nCon            
            legend(conNames)
        end
        set(bHan,'FaceColor',curvCatColors(j,:));
        set(bHan,'LineStyle',conStyles{k})        
        set(bHan,'LineWidth',3);
    end        
    
end
for j = 1:nCurvCat
    for k = 1:nCon
        currX = xEq(j,k);
        eHan = errorbar(currX,meanPer(k,j),semPer(k,j),'LineWidth',3);                
    end
    
    if hMW(k)
        starStr = '*';        
    else
        starStr = '';
    end    
    text(xEq(j,1),max(meanPer(:,j))+.1,[starStr ' p=' num2str(pMW(j),2)],'FontSize',14)        
end

set(currAx,'XTick',1:nCurvCat);
set(currAx,axProps{:})
xlim([.5 4.5])
ylim([0 1])

xlabel('Category #',labProps{:})
ylabel('Fraction of Cell Surface',labProps{:})

title({'Curvature category comparison',['WT n=' num2str(nMovPer(1)) ' Bleb n=' num2str(nMovPer(2))],'P Values from Mann-Whitney rank sum test','Error bars are +/- SEM'})

figName=  [outDir filesep 'curvature category histogram comparison WT and Bleb'];
if saveFigs,mfFigureExport(currFig,figName,figExpProps{:});end


%% ------------- Branch Varation around a mean figures ------ %%

%skPPFile = load('L:\nih\Post Processing Myo Inhib and WT\Skeleton Post Proc\WT\pruned skeleton post processing.mat');
skPPFile = load('W:\Hunter\orchestra_files_and_backup_merged\nih\Post Processing Myo Inhib and WT\Skeleton Post Proc\WT\pruned skeleton post processing.mat');

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

%parDir = 'L:\nih\Post Processing Myo Inhib and WT\Mask Geometry';
parDir = 'W:\Hunter\orchestra_files_and_backup_merged\nih\Post Processing Myo Inhib and WT\Mask Geometry';
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


%% ---------------- Branch angle vs. curvature   -------------- %%

% TOO BUMPY movPaths = {'W:\Hunter\orchestra_files_and_backup_merged\nih\Hunter data 2012_09_fixedcells\phall488-mem\488_dashPhal\488phal_dash_2\movieData.mat'};
% TOO BUMPY              W:\Hunter\orchestra_files_and_backup_merged\nih\act-mem_2013_01\act14\movieData.mat
movPaths = {'W:\Hunter\orchestra_files_and_backup_merged\nih\4D gfpMIIB fix and stain\control\control10\movieData.mat',...
    'W:\Hunter\orchestra_files_and_backup_merged\nih\4D gfpMIIB fix and stain\control\control12\movieData.mat',...
    'W:\Hunter\orchestra_files_and_backup_merged\nih\4D gfpMIIB fix and stain\control\control12\movieData.mat',...
    'W:\Hunter\orchestra_files_and_backup_merged\nih\4D gfpMIIB fix and stain\control\control13\movieData.mat',...
    'W:\Hunter\orchestra_files_and_backup_merged\nih\act-mem_2013_01\act11\movieData.mat',...
    'W:\Hunter\orchestra_files_and_backup_merged\nih\act-mem_2013_01\act15\movieData.mat',...
    'W:\Hunter\orchestra_files_and_backup_merged\nih\act-mem_2013_01\act17\movieData.mat',...
    'W:\Hunter\orchestra_files_and_backup_merged\nih\myoIIA60X_2103_02\set12\movieData.mat',...
    'W:\Hunter\orchestra_files_and_backup_merged\nih\myoIIA60X_2103_02\set13\movieData.mat'};


nMov = numel(movPaths);

useMicrons = true;


%%
curvRange = [.5 1.5];


closeCamPos = [-707.568 -1187.57 -3353.38;
               31.9917 -585.096 -2824.59;
               354.055 225.226 -2939.69;
               452.168 633.477 -3601.47;
               324.636 -38.1048 -3085.53;
               690 -357.371 3266.87;
               1021.91 -360.472 2240.19;
               1637.67 422.224 4450.42;
               315.364 -23.4243 -4214.31];
                  
closeCamAng = [1.23545 0.813574 0.738016 0.597469 1.19621 1.3278 0.911002 0.838506 .644051];
                

closeCamTarg = [291.317 321.576 50.3184;
                291.859 191.566 31.5454;
                354.055 225.226 31.5454;
                454.386 379.359 32.722;
                230.02 160.26 57.3853;
                301.913 207.299 43.3853;
                335.23 205.579 37.5954;
                383.93 301.503 57.8853;
                448.818 244.242 62.8853];
                
lPos = [-1 -1 -1;
        -1 -1 -1;
        -1 1 -1;
        -1 1 -1;
        -1 -1 -1;
        -1 -1 1;
        -1 -1 1;
        -1 -1 1;
         1  0 -1];
     
iFramePer = [1 1 1 1 1 1 1 1 1 1];
    

for iMov = 1:nMov

    MD = MovieData.load(movPaths{iMov},0);
    
    %TEMP
    %MD = ML.movies_{j}

    iMGProc =MD.getProcessIndex('MaskGeometry3DProcess',1,0);
    currMG = MD.processes_{iMGProc}.loadChannelOutput(1,iFramePer(iMov));
   
    currFig = fsFigure(.6); 
    pHan = showMaskSurfaceProp(currMG,'LAcurv');
   
    set(pHan,'FaceVertexCData',get(pHan,'FaceVertexCData') / MD.pixelSize_ * 1e3)
    
    caxis(curvRange)
    
    lh2 = light;
    set(lh2,'Position',lPos(iMov,:))
    %set(lh2,'Position',[1 0 -1])
    
    axHan = gca;
    
    set(axHan,'CameraPosition',closeCamPos(iMov,:))
    set(axHan,'CameraViewAngle',closeCamAng(iMov))
    set(axHan,'CameraTarget',closeCamTarg(iMov,:));

    figName = ['Branch angle and curvature example ' num2str(iMov)];
    if saveFigs,mfFigureExport(currFig,[outDir filesep figName]);end



end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% -------------- Branch # Vs. Curvature Figure Panel ------------ %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Get branch number and curvature stats for all movies
skPP = load('W:\Hunter\orchestra_files_and_backup_merged\nih\Post Processing Myo Inhib and WT\Skeleton Post Proc\All\pruned skeleton post processing.mat');
mgPP = load('W:\Hunter\orchestra_files_and_backup_merged\nih\Post Processing Myo Inhib and WT\Mask Geometry\All\combined stats and hists all movies.mat');
%And the condition and movie name information
load('W:\Hunter\orchestra_files_and_backup_merged\nih\movieListSuccesfullyProcessed and indices.mat');



%%
branchField = 'nTipsPerFrame';
branchName = '# Branch Tips';
curvField = 'mapcMedPerMov';
curvName = 'Per-Timepoint Median of Max Absolute Curvature Component';
curvConv = 1e3;%Flip sign and convert to microns
curvUnits = '1/microns';

%Get the movie index corresponding to each datapoint
nFramesPer = cellfun(@numel,skPP.(branchField));
nMov = numel(nFramesPer);
movIndex = arrayfun(@(x)(ones(1,nFramesPer(x))*x),1:nMov,'Unif',0);
movIndex = [movIndex{:}];

movCols = randomColormap(nMov,42);

x = [mgPP.(curvField){:}] * curvConv;
y = [skPP.(branchField){:}];

%% --- Plot including outliers for posterity --- %

currFig = figure;
scatter(x,y,5,movIndex,'filled')
colormap(movCols)
xlabel([curvName ', ' curvUnits])
ylabel(branchName)
title({'Curvature vs. branch number, including outliers',...
        'Color indicates movie, one spot per frame',...
        ['n = ' num2str(nMov) ' cells with m = ' num2str(sum(nFramesPer)) ' time points ']});

figName = [outDir filesep 'curvature vs branch number panel with outliers included'];
if saveFigs,mfFigureExport(currFig,figName,figExpProps{:});end

%% --- Plot excluding outliers colored by cell --- %

%Find and remove the outliers. Not that this is not completely
%statistically motivated - The first 62 are from two of the oldest movies
%which are consistent outliers in a variety of measures, and have a
%questionable (almost certainly incorrect, and in any case unverifiable)
%pixel size. The last 6 are from two where the segmented surface gets hairy
%(very noisy) very quickly due to bleaching and the frame selection didn't
%completely remove the artifact.
[iOutlier,iInlier] = detectOutliers([mgPP.mapcMedPerMov{:}],3);

x(iOutlier) = [];
y(iOutlier) = [];
movIndex(iOutlier) = [];
nPoints = numel(x);
nMov = numel(unique(movIndex));


currFig = figure;
scatter(x,y,5,movIndex,'filled')
colormap(movCols)
xlabel([curvName ', ' curvUnits])
ylabel(branchName)
title({'Curvature vs. branch number, outliers excluded',...
        'Color indicates movie, one spot per frame',...
        ['n = ' num2str(nMov) ' cells with m = ' num2str(sum(nFramesPer)) ' time points ']});

figName = [outDir filesep 'curvature vs branch number panel colored by cell'];
if saveFigs,mfFigureExport(currFig,figName,figExpProps{:});end

%% ----- Correlation Test ---- %%

%Test non-randomness of (nonlinear) correlation

[tauKendall,pValKendall] = corr(x',y','type','Kendall');
[rhoSpearman,pValSpearman] = corr(x',y','type','Spearman');

testString = ['Kendall''s Tau = ' num2str(tauKendall) ,' (p = ' num2str(pValKendall) '), Spearmann''s Rho = ' num2str(rhoSpearman) ' (p = ' num2str(pValSpearman) ')'];


%Since bob decided to take out the other conditions we re-calc the
%correlation with just the WT and blebb:
xWTBleb = x;
xWTBleb = xWTBleb(isWT(movIndex) | isBleb(movIndex));
yWTBleb = y;
yWTBleb = yWTBleb(isWT(movIndex) | isBleb(movIndex));
nPointsWTBleb = numel(xWTBleb);

[tauKendallWTBleb,pValKendallWTBleb] = corr(xWTBleb',yWTBleb','type','Kendall');
[rhoSpearmanWTBleb,pValSpearmanWTBleb] = corr(xWTBleb',yWTBleb','type','Spearman');
testStringWTBleb = ['Kendall''s Tau = ' num2str(tauKendallWTBleb) ,' (p = ' num2str(pValKendallWTBleb) '), Spearmann''s Rho = ' num2str(rhoSpearmanWTBleb) ' (p = ' num2str(pValSpearmanWTBleb) ')'];

%And for just the WT alone for good measure
xWT = x;
xWT = xWT(isWT(movIndex));
yWT = y;
yWT = yWT(isWT(movIndex));
nPointsWT = numel(xWT);

[tauKendallWT,pValKendallWT] = corr(xWT',yWT','type','Kendall');
[rhoSpearmanWT,pValSpearmanWT] = corr(xWT',yWT','type','Spearman');
testStringWT = ['Kendall''s Tau = ' num2str(tauKendallWT) ,' (p = ' num2str(pValKendallWT) '), Spearmann''s Rho = ' num2str(rhoSpearmanWT) ' (p = ' num2str(pValSpearmanWT) ')'];



%% -----Plot and corr colored by condition excluding outliers --- %


condNames = {'Untreated','Blebbistatin','C3','Y'};
condInd = zeros(nPoints,1);
condInd(isWT(movIndex)) = 1;
condInd(isBleb(movIndex)) = 2;
condInd(isC3(movIndex)) = 3;
condInd(isY(movIndex)) = 4;

condCols = [0 0 1;
            1 0 0;
            0 1 0;
            0 1 1];

currFig = figure;
hold on
for j = 1:size(condCols,1)
    plot(x(condInd == j),y(condInd == j),'.','color',condCols(j,:));
end
legend(condNames)
xlabel([curvName ', ' curvUnits])
ylabel(branchName)
title({'Curvature vs. branch number, outliers excluded , one spot per frame',...
        'Color indicates treatment',...
        ['n = ' num2str(nMov) ' cells with m = ' num2str(sum(nFramesPer)) ' time points '],testString});

figName = [outDir filesep 'curvature vs branch number panel colored by treatment'];
if saveFigs,mfFigureExport(currFig,figName,figExpProps{:});end

%% -----Plot and corr colored by condition excluding outliers Bleb and WT only --- %


condNames = {'Untreated','Blebbistatin'};
condInd = zeros(nPointsWTBleb,1);
condInd(isWT(movIndex)) = 1;
condInd(isBleb(movIndex)) = 2;

condCols = [0 0 1;
            1 0 0];            
currFig = figure;
hold on
for j = 1:size(condCols,1)
    plot(x(condInd == j),y(condInd == j),'.','color',condCols(j,:));
end

legend(condNames)
xlabel([curvName ', ' curvUnits])
ylabel(branchName)
title({'Curvature vs. branch number, outliers excluded , one spot per frame',...
        'Color indicates treatment',...
        ['n = ' num2str(nnz(isWT | isBleb)) ' cells with m = ' num2str(sum(nFramesPer(isWT | isBleb))) ' time points '],testStringWTBleb});

figName = [outDir filesep 'curvature vs branch number panel colored by treatment bleb and WT only'];
if saveFigs,mfFigureExport(currFig,figName,figExpProps{:});end

%% -----Plot and corr WT only --- %


    
currFig = figure;
hold on
plot(xWT,yWT,'b.')

legend('Untreated')
xlabel([curvName ', ' curvUnits])
ylabel(branchName)
title({'Curvature vs. branch number, outliers excluded , one spot per frame',...
        'Untreated only',...
        ['n = ' num2str(nnz(isWT)) ' cells with m = ' num2str(sum(nFramesPer(isWT))) ' time points '],testStringWT});

figName = [outDir filesep 'curvature vs branch number panel colored by treatment WT only'];
if saveFigs,mfFigureExport(currFig,figName,figExpProps{:});end

%% -----Plot and corr colored by treated/untreated only excluding outliers --- %

%To simplify things we just color by WT or myosin inhibited.
condNames = {'Untreated','Myosin Inhibited'};
condInd = ones(nPoints,1)*2;
condInd(isWT(movIndex)) = 1;

condCols = [0 0 1;
            1 0 0];

currFig = figure;
scatter(x,y,5,condInd,'filled')
colormap(condCols)
xlabel([curvName ', ' curvUnits])
ylabel(branchName)
title({'Curvature vs. branch number, outliers excluded , one spot per frame',...
        'Color indicates treated (blue untreated, red myosin inhibited)',...
        ['n = ' num2str(nMov) ' cells with m = ' num2str(sum(nFramesPer)) ' time points '],testString});

figName = [outDir filesep 'curvature vs branch number panel colored by treated or untreated'];
if saveFigs,mfFigureExport(currFig,figName,figExpProps{:});end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ------------ Curv Vs. Intensity Combined Corr Panel ----------- %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figDirs = {'W:\Hunter\orchestra_files_and_backup_merged\nih\myoIIA60X_2103_02\PostProcessing\Whole-Images First Frame Only\',...
           'W:\Hunter\orchestra_files_and_backup_merged\nih\4D gfpMIIB fix and stain\Intensity Analysis\WT\',...
            'W:\Hunter\orchestra_files_and_backup_merged\nih\act-mem_2013_01\Post Processing\Masked Intensity\'};
condNames = {'MIIA','MIIB','Actin'};
nCond = numel(condNames);

%figNameEnd = ' versus Mean Intensity channel 1 average of per-cell curv correlations plot by percentile.fig';
figNameEnd = ' versus Mean Intensity channel 1 plot by percentile.fig';

curvTypes = {'Mean Curvature','Gaussian Curvature','Max Absolute Curvature'};
nCurvTypes = numel(curvTypes);

condCols = [0 0 1 ;
            0 1 0 ;
            1 0 0 ];
x = cell(nCond,1);
y = cell(nCond,1);
ciH = cell(nCond,1);
ciL = cell(nCond,1);

for iCurv = 1:nCurvTypes
    
    currFig = figure;
    hold on
    
    for iCond = 1:nCond
    
        %Didn't save the combined data to .mat so we just extract it from
        %the figures
        figFile = [curvTypes{iCurv} figNameEnd];
        
        figHan = open([figDirs{iCond} filesep figFile]);
        axHan = get(figHan,'CurrentAxes');
        datHan = get(axHan,'Children');
        iChild = 4;%The first thing plotted is the last child, which is the data itself
        x{iCond} = get(datHan(iChild),'XData');
        y{iCond} = get(datHan(iChild),'YData');
        ciH{iCond} = get(datHan(2),'YData');
        ciL{iCond} = get(datHan(3),'YData');
        close(figHan);
        
        figure(currFig)
        plot(x{iCond},y{iCond},'color',condCols(iCond,:))
        xlabel('Fluorescence Intensity Percentile')
        ylabel([curvTypes{iCurv} ', 1/microns'])
       
        
    end
    legend(condNames)
    %Do CIs last so legend works
    for iCond = 1:nCond
        patch([x{iCond} x{iCond}(end:-1:1)],[ciH{iCond} ciL{iCond}(end:-1:1)],condCols(iCond,:),'EdgeColor','none','FaceAlpha',.2)
        
    end
    
    figName = [curvTypes{iCurv} ' vs  intensity all labels overlay per sample CI'];
    if saveFigs,mfFigureExport(currFig,[outDir filesep figName]);end     
    
end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% -------------- Spider Curvature Figure   ---------------------- %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



MD = MovieData.load('W:\Hunter\orchestra_files_and_backup_merged\nih\Spiders\Cell7c561_2x2\movieData.mat',0);
iFrame = 1;
intAn = MD.processes_{MD.getProcessIndex('MaskedIntensity3DProcess',1,0)}.loadChannelOutput(1);
maskGeo = MD.processes_{MD.getProcessIndex('MaskGeometry3DProcess',1,0)}.loadChannelOutput(1,iFrame);

%%


useMicrons = true;
%[curvTypes,curvNames,curvUnits,curvConv] = getCurveTypeFields(MD.pixelSize_,useMicrons);
%[intTypes,intNames] = getIntTypeFields;
curvTypes = {'LAmean','LAcurv'};
curvNames = {'Mean Curvature','Maximum Absolute Curvature Component'};
curvUnits = {'1/microns','1/microns'};
curvConv = 1e3/MD.pixelSize_ * [-1 1];

iCurvShow = [2 5];
cAxisPer = [-.2 0.2;
            0.33 2];
nCurvShow = numel(iCurvShow);
iIntShow = 1;

camPos = [-194.855 136.331 2006.79;
          -240.738 38.7671 -1899.49];
camTarg = [289.494 220.799 44.3678;
           284.764 218.61 45.8479];
camAng =[ 5.84887 5.59519];
nPos = numel(camAng);
viewName = {'Top view','bottom view'};

lPos = [-1 -1 1;
        -1 -1 -1];
for j = 1%:nCurvShow
    
    currFig = fsFigure(.75);
    
    pHan = showMaskSurfaceProp(maskGeo,curvTypes{j});
    set(pHan,'FaceVertexCData',get(pHan,'FaceVertexCData') * curvConv(j))
    colormap(jet(1024))
    caxis(cAxisPer(j,:))
    box off
    axis off
    axis ij    
    set(currFig,'color','w')
    lh2 = light;
    
    
    for k= 1%:nPos
        set(gca,'CameraPosition',camPos(k,:));
        set(gca,'CameraTarget',camTarg(k,:));
        set(gca,'CameraViewAngle',camAng(k));
        set(lh2,'Position',lPos(k,:))
        
        figName = ['Spider figure ' curvNames{j} ' ' viewName{k} 'jet'];
        if saveFigs,mfFigureExport(currFig,[outDir filesep figName]);end        
        
        colormap(curvColormap(256))
        
        figName = ['Spider figure ' curvNames{j} ' ' viewName{k} 'curvmap'];
        if saveFigs,mfFigureExport(currFig,[outDir filesep figName]);end        
        
        
    end
    
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% -------------- Example MII ROI Time Series -------------------- %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%TEMP- for sending to bob
outDir = 'W:\Hunter\orchestra_files_and_backup_merged\nih\figs for bob\9_3_2013\MII ROI Examples';


load('W:\Hunter\orchestra_files_and_backup_merged\nih\myoIIA60X_2103_02\movieListAllROIs and indices for processing.mat')
ML = MovieList.load('W:\Hunter\orchestra_files_and_backup_merged\nih\myoIIA60X_2103_02\movieListAllROIs.mat',0);

%%

iShow = find(isGood)';
nShow = numel(iShow);
iCurvShow = 2;
curvConv = 1e3;
curvUnits = '1/microns';
iIntShow = 1;

typeInd = ones(numel(isGood),1)*3;%We label protrusions by elimination
typeInd(iGoodFluct) = 1;
typeInd(iGoodRetract) = 2;
typeInd(~isGood) = 0;
typeNames = {'Fluctuation','Retraction','Protrusion'};

for j = iShow

    intAn = ML.movies_{j}.processes_{ML.movies_{j}.getProcessIndex('MaskedIntensity3DProcess',1,0)}.loadChannelOutput(1);

    tData = 0:ML.movies_{j}.timeInterval_:ML.movies_{j}.timeInterval_*(ML.movies_{j}.nFrames_-1);
    currFig = figure;
    [ah,h1,h2] = plotyy(tData,intAn.intStats.(['PerFrameMean' intAn.intTypes{iIntShow}]),...
                        tData,intAn.curvStats.(['PerFrameMean' intAn.curvTypes{iCurvShow}]) .* curvConv);
    set(get(ah(1),'Ylabel'),'String',['Photobleach Corrected ' intAn.intNames{iIntShow}])
    set(get(ah(2),'Ylabel'),'String',[intAn.curvNames{iCurvShow} ', ' curvUnits])
    xlabel('Time, Seconds')
    
    figName = ['Example '  typeNames{typeInd(j)} ' time series from ROI ' num2str(j) ];
    if saveFigs,mfFigureExport(currFig,[outDir filesep figName]);end        
    
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%% SUPPLMEMENTAL FIGURES %%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


outDir = 'W:\Hunter\orchestra_files_and_backup_merged\nih\Figures\Supplemental Figure Panels';
saveFigs = false;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ------------------------ Filtering ---------------------------- %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


movPath = 'W:\Hunter\orchestra_files_and_backup_merged\nih\data_2_2011_fast_and_c3\Faster\Faster3\Faster3a3\movieData.mat';
MD = MovieData.load(movPath,0);

iFrame = 1;
iChan = 1;
imName = MD.channels_(iChan).getImageFileNames(iFrame);
imName = imName{1};

%Get the raw images and filter response
im = stackRead([MD.channels_(iChan).channelPath_ filesep imName]);

%%


[maxResp,d2X,d2Y,d2Z] = multiscaleSurfaceFilter3D(im);

allImShow = cat(5,double(im),maxResp,d2X,d2Y,d2Z);
clear im;clear maxResp;clear d2X;clear d2Y,clear d2Z

%%



nImShow = size(allImShow,5);
imNames = {'Fluorescence','Filter Response Magnitude','Filter X Response','Filter Y Response','Filter Z Response'};
imCols = [1 0 0 ;
          1 0 0 ;
          0 1 0 ;
          0 0 1 ;
          1 1 0 ];
alpha = .5;

%%

camQuat = [0.2064   -0.2159   -0.6714    0.6782];%Orientation qauternion
camPos = [-13.6024   68.2360  101.5930];
camFoc = [131.7622];
camHt = [.7854];

nView = numel(camHt);

%%


[imsApp,ICEconn] = imarisShowArray(allImShow,[],[MD.pixelSize_ MD.pixelSize_ MD.zSpacing_]/1e3);


%%

pctileUse = .1;

useMax = false;
useBack = true;
for iIm = 1:nImShow
    currIm = allImShow(:,:,:,1,iIm);    
    [backMean,backSTD] = robustMean(double(currIm(:)),[],2);        
    if useMax        
        showRange(2) = max(currIm(:));        
    else
        showRange(2) = prctile(currIm(:),[100-pctileUse]);
    end    
    if useBack
       showRange(1) = backMean; 
    end        
    imsApp.GetDataSet.SetChannelRange(iIm-1,showRange(1),showRange(2));
    for k = 1:nView
        imsApp.GetSurpassCamera.SetFocus(camFoc(k));
        imsApp.GetSurpassCamera.SetHeight(camHt(k));
        imsApp.GetSurpassCamera.SetOrientationQuaternion(camQuat(k,:));
        imsApp.GetSurpassCamera.SetPosition(camPos(k,:));                            
        imsApp.GetDataSet.SetChannelColorRGBA(iIm-1,ICEconn.mapRgbaVectorToScalar([imCols(iIm,:),alpha]));        
    end 
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ------------------------ Segmentation ------------------------- %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%we reproduce everything here so we can show intermediate steps

%%
nSTDintensity = 3;
nSTDsurface = 2;
currIm = allImShow(:,:,:,1,1);
maxResp = allImShow(:,:,:,1,2);

[backMean,backSTD] = robustMean(double(currIm(:)),[],2);

backThresh = backMean + nSTDintensity*backSTD;
backMask = currIm > backThresh;

foreThresh = thresholdOtsu(currIm(backMask(:)));
foreMask = currIm>foreThresh;
surfBackMean = mean(maxResp(:));
surfBackSTD = std(maxResp(:));
surfThresh = surfBackMean + (nSTDsurface*surfBackSTD);
surfMask = maxResp > surfThresh;
mask = surfMask | foreMask;

ppp.MinVolume = 25;
ppp.NumObjects = 1; 
ppp.ClosureRadius = 2;
ppp.FillHoles = 2;
ppMask = postProcess3DMask(mask,ppp);

tmpMask = surfMask;
tmpMask(foreMask) = false;%To prevent Imaris from blending the colors;
allMaskShow = cat(5,single(currIm),single(backMask),single(foreMask),single(tmpMask),single(mask),single(ppMask));
%%

nImShow = size(allMaskShow,5);
[imsApp,ICEconn] = imarisShowArray(allMaskShow,[],[MD.pixelSize_ MD.pixelSize_ MD.zSpacing_]/1e3);

%%



imCols = [1 0 0 ;
          .5 .5 .5 ;
          0 0 1 ;
          1 1 0 ;
          1 0 0 ;
          0 1 0 ];
alpha = .5;


pctileUse = .1;

useMax = false;
useBack = true;
for iIm = 1:nImShow
    
    if iIm == 1
        currIm = allMaskShow(:,:,:,1,iIm);    
        [backMean,backSTD] = robustMean(double(currIm(:)),[],2);        
        if useMax        
            showRange(2) = max(currIm(:));        
        else
            showRange(2) = prctile(currIm(:),[100-pctileUse]);
        end    
        if useBack
           showRange(1) = backMean; 
        end        
    else
        showRange = [0 2];
    end
    imsApp.GetDataSet.SetChannelRange(iIm-1,showRange(1),showRange(2));
    for k = 1:nView
        imsApp.GetSurpassCamera.SetFocus(camFoc(k));
        imsApp.GetSurpassCamera.SetHeight(camHt(k));
        imsApp.GetSurpassCamera.SetOrientationQuaternion(camQuat(k,:));
        imsApp.GetSurpassCamera.SetPosition(camPos(k,:));                            
        imsApp.GetDataSet.SetChannelColorRGBA(iIm-1,ICEconn.mapRgbaVectorToScalar([imCols(iIm,:),alpha]));        
    end 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ---------------- Surface creation / curv calc ---------------- %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

outDir = 'W:\Hunter\orchestra_files_and_backup_merged\nih\Figures\Supplemental Figure Panels\Surface Mesh and Curvature\';

%again we reproduce everything here to show intermediate steps
mask = make3DImageVoxelsSymmetric(ppMask,MD.pixelSize_,MD.zSpacing_);

maskFilt = fastGauss3D(mask,1,[5 5 5],2);

surf = isosurface(maskFilt,.25);

norms = isonormals(maskFilt,surf.vertices);        

[K,H] = surfaceCurvature(surf,norms);

H = -H;%for the new sign convention;

k1 = H + sqrt(H .^2 - K);
k2 = H - sqrt(H .^2 - K);

mp.GaussianCurvature = K;
mp.MeanCurvature = H;
mp.SmoothedSurface = surf;
mp.SurfaceNormals = norms;
mp.CurvaturePC1 = k1;
mp.CurvaturePC2 = k2;


locAvgCurv = calcLocalAvgCurvatures(mp,2e3/MD.pixelSize_,1,1);

distX = bwdist(~mask);

%% Matlab figures for surface gen and curv calc 




camPos = [653.083 -1246.1 1463.44;
          654.12 -1222.82 1463.7];
camTarg = [163.322 261.235 36.3869;
           164.359 284.52 36.6453];
camAng = [1.9496 .243723];
nView = numel(camAng);

%% Triangulated mesh & normals illustration
currFig = fsFigure(.75);
pHan = patch(mp.SmoothedSurface,'FaceColor',[.8 .8 .8],'FaceAlpha',1,...
    'EdgeColor','k');
axis ij
axis image
axis off
set(currFig,'color','w')
hold on
qHan = quiver3(mp.SmoothedSurface.vertices(:,1),mp.SmoothedSurface.vertices(:,2),mp.SmoothedSurface.vertices(:,3),...
              mp.SurfaceNormals(:,1),mp.SurfaceNormals(:,2),mp.SurfaceNormals(:,3),.75);
    
lighting flat
light
lh2 = light;
set(lh2,'Position',[0 1 -1])


for iView = 1:nView
    
    
    set(gca,'CameraPosition',camPos(iView,:));
    set(gca,'CameraTarget',camTarg(iView,:));
    set(gca,'CameraViewAngle',camAng(iView));



    figName = ['Surface mesh and normals illustration figure ' num2str(iView)];
    if saveFigs,mfFigureExport(currFig,[outDir filesep figName]);end        
    
end

%% Curvature illustration

currFig = fsFigure(.75);
pHan = patch(mp.SmoothedSurface,'FaceColor','flat','FaceAlpha',1,...
    'EdgeColor','k');
axis ij
axis image
axis off
set(currFig,'color','w')
hold on
% qHan = quiver3(mp.SmoothedSurface.vertices(:,1),mp.SmoothedSurface.vertices(:,2),mp.SmoothedSurface.vertices(:,3),...
%               mp.SurfaceNormals(:,1),mp.SurfaceNormals(:,2),mp.SurfaceNormals(:,3),.75);
    
lighting flat
light
lh2 = light;
set(lh2,'Position',[0 1 -1])

curvDat = [max(abs([k1 k2]),[],2) locAvgCurv.LocMeanMaxAbsCurvature];
curvNames = {'Raw Max','Locally Averaged Max'};
nCurv = size(curvDat,2);

pSat = 5;

for iCurv = 1:nCurv 
    
    cRange = prctile(curvDat(:,iCurv),[pSat 100-pSat]);
    
    for iView = 1:nView

        set(pHan,'FaceVertexCData',curvDat(:,iCurv));

        set(gca,'CameraPosition',camPos(iView,:));
        set(gca,'CameraTarget',camTarg(iView,:));
        set(gca,'CameraViewAngle',camAng(iView));

        caxis(cRange)

        figName = ['Surface curvature ' curvNames{iCurv} ' illustration figure ' num2str(iView)];
        if saveFigs,mfFigureExport(currFig,[outDir filesep figName]);end        

    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ------------ Surface & int corr / depth norm ------------------ %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

outDir = 'W:\Hunter\orchestra_files_and_backup_merged\nih\Figures\Supplemental Figure Panels\Curvature Correlation and Depth Norm';

%Use an image with myosin and actin 

%Has good cortical act and myo but membrane stain poor
movPath = 'W:\Hunter\orchestra_files_and_backup_merged\nih\Hunter data 2012_09_fixedcells\gfpmyoIIA\myoIIA_8\movieData.mat';


%movPath = 'W:\Hunter\orchestra_files_and_backup_merged\nih\4D gfpMIIB fix
%and stain\control\control10\movieData.mat';%NOT STRONGLY CORTICAL
MD = MovieData.load(movPath,0);



%%

iFrame = 1;
chanShow = 1:numel(MD.channels_);
nChanShow= numel(chanShow);
clear im
for iChan = 1:nChanShow
    imName = MD.channels_(chanShow(iChan)).getImageFileNames(iFrame);
    imName = imName{1};
    im(:,:,:,iChan) = stackRead([MD.channels_(chanShow(iChan)).channelPath_ filesep imName]);
end
imSize = size(im);

%%

maskName = MD.processes_{MD.getProcessIndex('SegmentationProcess3D',1,0)}.getOutMaskFileNames(1);%Only one channel of masks, and they are associated with chan 1
maskPath = MD.processes_{MD.getProcessIndex('SegmentationProcess3D',1,0)}.outFilePaths_{1};
maskName = maskName{1};

m = tif3Dread([maskPath filesep maskName{iFrame}]);

distX = bwdist(~m);
iChanDN = 1;
imDN = depthNormalizeImage(im(:,:,:,iChanDN),distX,false(imSize(1:3)));

%% ---- Example Raw and Depth Normalized Images ----- %%



iSlice= 16;
xl = [290 615];
yl = [405 648];
satPct = 1;
scBarSize = 5;%Size in microns
nPixScale = scBarSize*1e3/MD.pixelSize_;
scaleStr = [num2str(scBarSize) ' microns'];

% ---- Example Image before depth normalization ---- %

currFig = fsFigure(.6);
imHan = imshow(im(:,:,iSlice),[]);
xlim(xl);
ylim(yl);
plotScaleBar(nPixScale,'Label',scaleStr,'Location','SouthEast')
figName = 'Depth normalization example MIIA-GFP image before normalization';
if saveFigs,mfFigureExport(currFig,[outDir filesep figName]);end        

%saturateImageColormap(gca,satPct);

% ---- Example Image after depth normalization ---- %

currFig = fsFigure(.6);
imHan = imshow(imDN(:,:,iSlice),[]);
xlim(xl);
ylim(yl);
plotScaleBar(nPixScale,'Label',scaleStr,'Location','SouthEast')
figName = 'Depth normalization example MIIA-GFP image after normalization';
if saveFigs,mfFigureExport(currFig,[outDir filesep figName]);end        



%%  --- Curvature comparisons with and without depth normalization  --- %%

%I didn't save this data to .mat file so just open the figures and extract the data to overlay

figDirs = {'W:\Hunter\orchestra_files_and_backup_merged\nih\myoIIA60X_2103_02\PostProcessing\Whole-Images First Frame Only\',...
           'W:\Hunter\orchestra_files_and_backup_merged\nih\4D gfpMIIB fix and stain\Intensity Analysis\WT\',...
            'W:\Hunter\orchestra_files_and_backup_merged\nih\act-mem_2013_01\Post Processing\Masked Intensity\'};

curvTypes = {'Mean Curvature','Gaussian Curvature','Max Absolute Curvature'};
nCurvTypes = numel(curvTypes);
condNames = {'MIIA','MIIB','Actin'};
nCond = numel(condNames);

figNameEnd = ' average of per-cell curv correlations plot by percentile.fig';
%figNameEnd = ' plot by percentile.fig';

for iCond = 1:numel(figDirs);

    for iCurv = 1:nCurvTypes
        % - Get the raw data without depth norm

        figFile = [curvTypes{iCurv} ' versus Mean Intensity channel 1' figNameEnd];
        figHan = open([figDirs{iCond} figFile]);
        axHan = get(figHan,'CurrentAxes');
        datHan = get(axHan,'Children');
        iChild = 4;%The first thing plotted is the last child, which is the data itself
        xRaw = get(datHan(iChild),'XData');
        yRaw = get(datHan(iChild),'YData');
        ciRawH = get(datHan(2),'YData');
        ciRawL = get(datHan(3),'YData');
        close(figHan);
        
        % - Get the depth normalized data

        figFile = [curvTypes{iCurv} ' versus Depth Normalized Mean Intensity channel 1' figNameEnd];
        figHan = open([figDirs{iCond} figFile]);
        axHan = get(figHan,'CurrentAxes');
        datHan = get(axHan,'Children');
        iChild = 4;%The first thing plotted is the last child, which is the data itself
        xDN = get(datHan(iChild),'XData');
        yDN = get(datHan(iChild),'YData');
        ciDNH = get(datHan(2),'YData');
        ciDNL = get(datHan(3),'YData');
        close(figHan);

        %%

        currFig = figure;
        hold on
        plot(xRaw,yRaw,'b')
        plot(xDN,yDN,'r')
        legend('Raw Fluorescence','Depth Normalized Fluorescence')
        xlabel([condNames{iCond} ' Fluorescence Intensity Percentile'])
        ylabel(curvTypes{iCurv})
        patch([xRaw xRaw(end:-1:1)],[ciRawH ciRawL(end:-1:1)],'b','EdgeColor','none','FaceAlpha',.2)
        patch([xDN xDN(end:-1:1)],[ciDNH ciDNL(end:-1:1)],'r','EdgeColor','none','FaceAlpha',.2)
        xlim([0.5 99.5])

        figName = [curvTypes{iCurv} ' vs ' condNames{iCond} ' intensity depth normalization overlay per cell CI'];
        if saveFigs,mfFigureExport(currFig,[outDir filesep figName]);end        
    end
end

%% ----- Artifactual correlations and depth norm illustration ---- %%

%We just need the geometry, so pick one with relatively small image so the
%calc doesn't take forever
movPath = 'W:\Hunter\orchestra_files_and_backup_merged\nih\4D gfpMIIB fix and stain\control\control10\movieData.mat';
%movPath = 'W:\Hunter\orchestra_files_and_backup_merged\nih\20090605_LatestData_Blood\TdCAAX2\set3\pos4\movieData.mat';

MD = MovieData.load(movPath,0);

%% get mask and setup test image

iFrame = 1;
iSlice = 48;%Slice for showing cross-sections
maskName = MD.processes_{MD.getProcessIndex('SegmentationProcess3D',1,0)}.getOutMaskFileNames(1);%Only one channel of masks, and they are associated with chan 1
maskPath = MD.processes_{MD.getProcessIndex('SegmentationProcess3D',1,0)}.outFilePaths_{1};
maskName = maskName{1};

m = tif3Dread([maskPath filesep maskName{iFrame}]);

m = make3DImageVoxelsSymmetric(m,MD.pixelSize_,MD.zSpacing_);


%% ---- create depth-only test image ----- %%


% Depth - only variation image 
testImDepth = bwdist(~m);
nSig = 1;
gSig = 4;
amp = 10;
testImDepth = amp*exp(-((testImDepth - gSig).^2) / gSig) + nSig*randn(size(testImDepth));
testImDepth(m) = testImDepth(m) + amp/2;%Add an offset in the cell interior to avoid dividing by small numbers in the depth normalization

%% Example image from depth only test image

currFig = figure;
tmp = testImDepth;
%For visualization we set the voxels outside the cell to zero (this doesn't
%matter for the correlation because they aren't being sampled anyways)
tmp(~m) = 0;
imshow(tmp(:,:,iSlice),[])
figName = ['example plane from depth-only varying test image nSig ' num2str(nSig) ' gSig ' num2str(gSig) ' amp ' num2str(amp)];

if saveFigs,mfFigureExport(currFig,[outDir filesep figName]);end        


%%  run depth norm and int vs. curv corr on the test image 

%testImDN = depthNormalizeImage(testImDepth,m,false(size(m)));

mp = MD.processes_{MD.getProcessIndex('MaskGeometry3DProcess',1,0)}.loadChannelOutput(1,iFrame);
sk = MD.processes_{MD.getProcessIndex('SkeletonPruningProcess',1,0)}.loadChannelOutput(1,iFrame);


bp = analyze3DImageMaskedIntensities(testImDepth,m,sk,mp,2e3/MD.pixelSize_,false(size(m)));
%bpRand = analyze3DImageMaskedIntensities(testImRand,m,sk,mp,2e3/MD.pixelSize_,false(size(m)));
%bpLat = analyze3DImageMaskedIntensities(testImLat,m,sk,mp,2e3/MD.pixelSize_,false(size(m)));

%%

intField = 'intForCurvSampMean';
intName = 'Intensity';
intFieldDN = 'intForCurvSampMeanDepthNorm';
intNameDN = 'Depth Normalized Intensity';

%curvField = 'meanCurvSampMean';
curvField = 'MaxAbsPCCurvSampMean';
curvName = 'Max Absolute Principle Curvature';
curvConv = 1e3/MD.pixelSize_;



alpha = .05;

maxVal= max(testImDepth(:));
nBins = 100;
ptilesUse = linspace(0,100,nBins);
histBins = prctile(bp.(intField)(:),ptilesUse);
histBinsDN = prctile(bp.(intFieldDN)(:),ptilesUse);
%histBins = linspace(min(bp.(intField)(:)),max(bp.(intField)(:)),nBins);
%histBinsBN = linspace(min(bp.(intFieldDN)(:)),max(bp.(intFieldDN)(:)),nBins);
meanCurv= nan(nBins-1,1);
curvCI = nan(nBins-1,2);
meanCurvDN= nan(nBins-1,1);
curvCIDN = nan(nBins-1,2);
for j = 1:(nBins-1)    
    currSamp = bp.(intField) > histBins(j) & bp.(intField) <= histBins(j+1);    
    meanCurv(j) = mean(bp.(curvField)(currSamp) * curvConv);        
    if nnz(currSamp)>1
        bootSamp = bootstrp(1e3,@mean,bp.(curvField)(currSamp) * curvConv);
        curvCI(j,:) = prctile(bootSamp,[alpha/2 100-alpha/2]);
    end
    
    currSampDN = bp.(intFieldDN) > histBinsDN(j) & bp.(intFieldDN) <= histBinsDN(j+1);    
    meanCurvDN(j) = mean(bp.(curvField)(currSampDN) * curvConv);        
    if nnz(currSampDN)>1
        bootSamp = bootstrp(1e3,@mean,bp.(curvField)(currSampDN) * curvConv);
        curvCIDN(j,:) = prctile(bootSamp,[alpha/2 100-alpha/2]);
    end
    
end
    
%% ---- Artifactual int vs. curv curves depth varying only --- %%

curvX = ptilesUse(1:end-1);
currFig = figure;
plot(curvX,meanCurv)
hold on
plot(curvX,meanCurvDN,'r')
xlabel('Intensity Percentile')
ylabel([curvName ', 1/microns'])
%%Doesn't work because it has NaNs
%patch([curvX curvX(end:-1:1)],[curvCI(:,1)',curvCI(end:-1:1,2)'],'b')
legend(intName,intNameDN);

plot(curvX,curvCI(:,1),'--b')
plot(curvX,curvCI(:,2),'--b')
plot(curvX,curvCIDN(:,1),'--r')
plot(curvX,curvCIDN(:,2),'--r')

xlim([min(curvX) max(curvX)])


figName = ['Artifactual correlation and depth norm ' curvName ' vs  int depth only varying nSig ' num2str(nSig) ' gSig ' num2str(gSig) ' amp ' num2str(amp)];
if saveFigs,mfFigureExport(currFig,[outDir filesep figName]);end        

%% ---- Artifactual int vs. curv 2D Hist --- %%

[N,C] = hist3([bp.(intFieldDN) bp.(curvField) * curvConv],[nBins nBins]);

currFig = figure;
%imagesc(C{1},C{2},log10(N))
imagesc(log10(N))


figName = ['Artifactual correlation and depth norm ' curvName ' vs  distance transform correlation noise sigma ' num2str(nSig)];
if saveFigs,mfFigureExport(currFig,[outDir filesep figName]);end        

%% lateral variation only image

latSig = 10*nSig;
outFile = [outDir filesep 'lateral variation test image lat sig ' num2str(latSig) '.mat'];

useSaved = false;
if ~useSaved
    
    
    %Get distance transform label matrix
    [distX,distL] = bwdist(~m);
    distL(~m) = 0;
    %And find all surface voxels
    surfVox = unique(distL(distX(:) == 1));
    nSurfVox = numel(surfVox);
    testImLat = zeros(size(m));
    randVals = latSig*randn(nSurfVox,1);
    tic
    for j = 1:nSurfVox
        testImLat(distL == surfVox(j)) = randVals(j);        
    end
    toc
    
    save(outFile,'testImLat','latSig','distX','distL','m','surfVox','nSurfVox','randVals');
else%This took like a half hour to generate so we save to disk
    
    load(outFile)
    
end

%% random image with and wihtout spatial autocorrelation

nRep = 10;
bpRandAC = cell(nRep,1);
doRand = false;
if doRand
    bpRand = cell(nRep,1);
end

for iRep = 1:nRep
    tic
    %Spatially autocorrelated random image
    testImRandAC = randn(size(m));
    acSig = 1;
    testImRandAC = filterGauss3D(testImRandAC,acSig);
    %And set mean and std equal to the IID random image
    testImRandAC = testImRandAC - mean(testImRandAC(:));
    testImRandAC = testImRandAC ./ std(testImRandAC(:)) * nSig;
    testImRandAC = testImRandAC + amp/2;
    
    bpRandAC{iRep} = analyze3DImageMaskedIntensities(testImRandAC,m,sk,mp,2e3/MD.pixelSize_,false(size(m)));
    
    if doRand
        %% No autocorrelation random image
        testImRand = nSig*randn(size(m));
        %testImRand(m) = testImRand(m) + amp/2;%Add an offset in the cell interior to avoid dividing by small numbers in the depth normalization
        testImRand = testImRand + amp/2;%Add an offset to the whole image for spatial autocorr calc

        bpRand{iRep} = analyze3DImageMaskedIntensities(testImRand,m,sk,mp,2e3/MD.pixelSize_,false(size(m)));
    end
    disp(iRep);
    toc
end
save([outDir filesep 'random and autocorr random test images and corr acSig ' num2str(acSig)],'bpRandAC','bpRand','acSig','amp','nSig')
bpRandAll = vertcat(bpRand{:});
bpRandACAll = vertcat(bpRandAC{:});

%% Example image from random test image

currFig = figure;
tmp = testImRand;
%For visualization we set the voxels outside the cell to zero (this doesn't
%matter for the correlation because they aren't being sampled anyways)
tmp(~m) = 0;
imshow(tmp(:,:,iSlice),[])
figName = ['example plane random test image nSig ' num2str(nSig)];

if saveFigs,mfFigureExport(currFig,[outDir filesep figName]);end        

%% Example image from autocorrelated random test image

currFig = figure;
tmp = testImRandAC;
%For visualization we set the voxels outside the cell to zero (this doesn't
%matter for the correlation because they aren't being sampled anyways)
tmp(~m) = 0;
imshow(tmp(:,:,iSlice),[])
figName = ['example plane autocorrelated random test image nSig ' num2str(nSig) ' ac sig ' num2str(acSig)];

if saveFigs,mfFigureExport(currFig,[outDir filesep figName]);end    


%% int vs. curv corr on random test image

intField = 'intForCurvSampMean';
intName = 'Intensity';
intFieldDN = 'intForCurvSampMeanDepthNorm';
intNameDN = 'Depth Normalized Intensity';

%curvField = 'meanCurvSampMean';
curvField = 'MaxAbsPCCurvSampMean';
curvName = 'Max Absolute Principle Curvature';
curvConv = 1e3/MD.pixelSize_;



alpha = .05;

nBins = 100;
ptilesUse = linspace(0,100,nBins);

histBins = prctile(vertcat(bpRandAll(:).(intField)),ptilesUse);
histBinsDN = prctile(vertcat(bpRandAll(:).(intFieldDN)),ptilesUse);

meanCurvRand= nan(nBins-1,1);
curvCIRand = nan(nBins-1,2);
meanCurvDNRand= nan(nBins-1,1);
curvCIDNRand = nan(nBins-1,2);
meanSampSizeRand = nan(nBins-1,1);
sampSizeCIRand = nan(nBins-1,2);

meanCurvRandAll= nan(nBins-1,nRep);
meanCurvDNRandAll= nan(nBins-1,nRep);
meanSampSizeRandAll = nan(nBins-1,nRep);


for j = 1:(nBins-1)    
    
    for k = 1:nRep
        currSamp = bpRandAll(k).(intField) > histBins(j) & bpRandAll(k).(intField) <= histBins(j+1);    
        meanCurvRandAll(j,k) = mean(bpRandAll(k).(curvField)(currSamp) * curvConv);        

        meanSampSizeRandAll(j,k) = mean(bpRandAll(k).nPixPerCurvSamp(currSamp));
        
        currSampDN = bpRandAll(k).(intFieldDN) > histBinsDN(j) & bpRandAll(k).(intFieldDN) <= histBinsDN(j+1);    
        meanCurvDNRandAll(j,k) = mean(bpRandAll(k).(curvField)(currSampDN) * curvConv);        
        
        
    end
    
    meanCurvRand(j) = mean(meanCurvRandAll(j,:));
    bootSamp =  bootstrp(1e3,@mean, meanCurvRandAll(j,:));
    curvCIRand(j,:) = prctile(bootSamp,[alpha/2 100-alpha/2]);
    
    meanSampSizeRand(j) = mean(meanSampSizeRandAll(j,:));
    bootSamp =  bootstrp(1e3,@mean, meanSampSizeRandAll(j,:));
    sampSizeCIRand(j,:) = prctile(bootSamp,[alpha/2 100-alpha/2]);
    
    meanCurvDNRand(j) = mean(meanCurvDNRandAll(j,:));
    bootSamp =  bootstrp(1e3,@mean, meanCurvDNRandAll(j,:));
    curvCIDNRand(j,:) = prctile(bootSamp,[alpha/2 100-alpha/2]);
    
end

%% ---- Artifactual int vs. curv curves random image --- %%

currFig = figure;
curvX = ptilesUse(1:end-1);
plot(curvX,meanCurvRand)
hold on
plot(curvX,meanCurvDNRand,'r')


xlabel('Intensity Percentile')
ylabel([curvName ', 1/microns'])
%%Doesn't work because it has NaNs
%patch([curvX curvX(end:-1:1)],[curvCI(:,1)',curvCI(end:-1:1,2)'],'b')
meanSampCurv = nanmean(vertcat(bpRandAll(:).(curvField)) * curvConv);
plot(xlim,[1 1]*meanSampCurv,'k')

legend(intName,intNameDN,'Mean');

plot(curvX,curvCIRand(:,1),'--b')
plot(curvX,curvCIRand(:,2),'--b')
plot(curvX,curvCIDNRand(:,1),'--r')
plot(curvX,curvCIDNRand(:,2),'--r')

xlim([min(curvX) max(curvX)])


figName = ['Artifactual correlation and depth norm ' curvName ' vs  random image noise sigma ' num2str(nSig)];
if saveFigs,mfFigureExport(currFig,[outDir filesep figName]);end        

%% ---- Sample size vs. intensity on random image ---- %%


currFig = figure;
curvX = ptilesUse(1:end-1);
plot(curvX,meanSampSizeRand)
hold on
hemiVol = 4/3 * pi * (2e3/MD.pixelSize_)^3 / 2;%Volume of hemisphere w/ sampling radius
%plot(xlim,[1 1]*hemiVol,'--k')
xlabel('Intensity Percentile')
ylabel('Mean Sample Size (voxels)')

figName = ['Sample size vs intensity percentile in random image'];
if saveFigs,mfFigureExport(currFig,[outDir filesep figName]);end        

%% ---- Sample size vs. curvature ---- %%


%This readout doesn't depend on the image intensities so just use one of
%the replicates

x = bpRand{1}.MaxAbsPCCurvSampMean;
y = bpRand{1}.nPixPerCurvSamp;

currFig = figure;
hold on
plot(x,y,'.')
[tau,pVal] = corr(x,y,'Type','Kendall');
xlabel('Maximum Absolute Curvature Component, 1/microns')
ylabel('Sample size, voxels')

figName = ['Sample size vs max abs curvature'];
if saveFigs,mfFigureExport(currFig,[outDir filesep figName]);end    

%% ---- Sample depth vs. curvature ---- %%


%This readout doesn't depend on the image intensities so just use one of
%the replicates

x = bpRand{1}.MaxAbsPCCurvSampMean;
y = bpRand{1}.meanDepthPerCurvSamp;

currFig = figure;
hold on
plot(x,y,'.')
[tau,pVal] = corr(x,y,'Type','Kendall');
xlabel('Maximum Absolute Curvature Component, 1/microns')
ylabel('Mean Sample Depth, voxels')

figName = ['Sample depth vs max abs curvature'];
if saveFigs,mfFigureExport(currFig,[outDir filesep figName]);end   

%% int vs. curv corr on autocorrelated random test image

intField = 'intForCurvSampMean';
intName = 'Intensity';
intFieldDN = 'intForCurvSampMeanDepthNorm';
intNameDN = 'Depth Normalized Intensity';

%curvField = 'meanCurvSampMean';
curvField = 'MaxAbsPCCurvSampMean';
curvName = 'Max Absolute Principle Curvature';
curvConv = 1e3/MD.pixelSize_;



alpha = .05;

maxVal= max(testImDepth(:));
nBins = 100;
ptilesUse = linspace(0,100,nBins);

histBins = prctile(vertcat(bpRandACAll(:).(intField)),ptilesUse);
histBinsDN = prctile(vertcat(bpRandACAll(:).(intFieldDN)),ptilesUse);

meanCurvRandAC= nan(nBins-1,1);
curvCIRandAC = nan(nBins-1,2);
meanCurvDNRandAC= nan(nBins-1,1);
curvCIDNRandAC = nan(nBins-1,2);

meanCurvRandAllAC= nan(nBins-1,nRep);
meanCurvDNRandAllAC= nan(nBins-1,nRep);



for j = 1:(nBins-1)    
    
    for k = 1:nRep
        currSamp = bpRandACAll(k).(intField) > histBins(j) & bpRandACAll(k).(intField) <= histBins(j+1);    
        meanCurvRandAllAC(j,k) = mean(bpRandACAll(k).(curvField)(currSamp) * curvConv);        

        currSampDN = bpRandACAll(k).(intFieldDN) > histBinsDN(j) & bpRandACAll(k).(intFieldDN) <= histBinsDN(j+1);    
        meanCurvDNRandAllAC(j,k) = mean(bpRandACAll(k).(curvField)(currSampDN) * curvConv);        
    end
    
    meanCurvRandAC(j) = mean(meanCurvRandAllAC(j,:));
    bootSamp =  bootstrp(1e3,@mean, meanCurvRandAllAC(j,:));
    curvCIRandAC(j,:) = prctile(bootSamp,[alpha/2 100-alpha/2]);
    
    meanCurvDNRandAC(j) = mean(meanCurvDNRandAllAC(j,:));
    bootSamp =  bootstrp(1e3,@mean, meanCurvDNRandAllAC(j,:));
    curvCIDNRandAC(j,:) = prctile(bootSamp,[alpha/2 100-alpha/2]);
    
end

%% ---- Artifactual int vs. curv curves random autocorr --- %%

currFig = figure;
curvX = ptilesUse(1:end-1);
plot(curvX,meanCurvRandAC)
hold on
plot(curvX,meanCurvDNRandAC,'r')
meanSampCurv = nanmean(vertcat(bpRandACAll(:).(curvField)) * curvConv);
plot(xlim,[1 1]*meanSampCurv,'k')

xlabel('Intensity Percentile')
ylabel([curvName ', 1/microns'])
%%Doesn't work because it has NaNs
%patch([curvX curvX(end:-1:1)],[curvCI(:,1)',curvCI(end:-1:1,2)'],'b')
legend(intName,intNameDN,'Mean');

plot(curvX,curvCIRandAC(:,1),'--b')
plot(curvX,curvCIRandAC(:,2),'--b')
plot(curvX,curvCIDNRandAC(:,1),'--r')
plot(curvX,curvCIDNRandAC(:,2),'--r')

xlim([min(curvX) max(curvX)])




figName = ['Artifactual correlation and depth norm ' curvName ' vs  autocorr random nSig ' num2str(nSig) ' ac sig ' num2str(acSig)];
if saveFigs,mfFigureExport(currFig,[outDir filesep figName]);end     


%% random image with spatial autocorrelation and noise

nRep = 10;
bpRandACN = cell(nRep,1);
eSig = nSig/2;%Sigma of IID noise added on top
acSig = 2;
for iRep = 1:nRep
    
    tic
    %Spatially autocorrelated random image
    testImRandACN = randn(size(m));    
    testImRandACN = filterGauss3D(testImRandACN,acSig);
    %And set mean and std equal to the IID random image
    testImRandACN = testImRandACN - mean(testImRandACN(:));
    testImRandACN = testImRandACN ./ std(testImRandACN(:)) * nSig;
    testImRandACN = testImRandACN + amp/2;
    testImRandACN = testImRandACN + randn(size(m)) * eSig;
    
    bpRandACN{iRep} = analyze3DImageMaskedIntensities(testImRandACN,m,sk,mp,2e3/MD.pixelSize_,false(size(m)));
            
    disp(iRep);
    toc
    
end

save([outDir filesep 'random autocorr with noise test images and corr acSig ' num2str(acSig) ' eSig ' num2str(eSig) '.mat'],'bpRandACN','acSig','amp','nSig','eSig')
bpRandACNAll = vertcat(bpRandACN{:});

%% int vs. curv corr on autocorrelated plus noise random test image

intField = 'intForCurvSampMean';
intName = 'Intensity';
intFieldDN = 'intForCurvSampMeanDepthNorm';
intNameDN = 'Depth Normalized Intensity';

%curvField = 'meanCurvSampMean';
curvField = 'MaxAbsPCCurvSampMean';
curvName = 'Max Absolute Principle Curvature';
curvConv = 1e3/MD.pixelSize_;



alpha = .05;
nBins = 100;
ptilesUse = linspace(0,100,nBins);

histBins = prctile(vertcat(bpRandACNAll(:).(intField)),ptilesUse);
histBinsDN = prctile(vertcat(bpRandACNAll(:).(intFieldDN)),ptilesUse);

meanCurvRandACN= nan(nBins-1,1);
curvCIRandACN = nan(nBins-1,2);
meanCurvDNRandACN= nan(nBins-1,1);
curvCIDNRandACN = nan(nBins-1,2);

meanCurvRandAllACN= nan(nBins-1,nRep);
meanCurvDNRandAllACN= nan(nBins-1,nRep);



for j = 1:(nBins-1)    
    
    for k = 1:nRep
        currSamp = bpRandACNAll(k).(intField) > histBins(j) & bpRandACNAll(k).(intField) <= histBins(j+1);    
        meanCurvRandAllACN(j,k) = mean(bpRandACNAll(k).(curvField)(currSamp) * curvConv);        

        currSampDN = bpRandACNAll(k).(intFieldDN) > histBinsDN(j) & bpRandACNAll(k).(intFieldDN) <= histBinsDN(j+1);    
        meanCurvDNRandAllACN(j,k) = mean(bpRandACNAll(k).(curvField)(currSampDN) * curvConv);        
    end
    
    meanCurvRandACN(j) = mean(meanCurvRandAllACN(j,:));
    bootSamp =  bootstrp(1e3,@mean, meanCurvRandAllACN(j,:));
    curvCIRandACN(j,:) = prctile(bootSamp,[alpha/2 100-alpha/2]);
    
    meanCurvDNRandACN(j) = mean(meanCurvDNRandAllACN(j,:));
    bootSamp =  bootstrp(1e3,@mean, meanCurvDNRandAllACN(j,:));
    curvCIDNRandACN(j,:) = prctile(bootSamp,[alpha/2 100-alpha/2]);
    
end

%% ---- Artifactual int vs. curv curves random autocorr plus noise --- %%

currFig = figure;
curvX = ptilesUse(1:end-1);
plot(curvX,meanCurvRandACN)
hold on
plot(curvX,meanCurvDNRandACN,'r')
meanSampCurv = nanmean(vertcat(bpRandACNAll(:).(curvField)) * curvConv);
plot(xlim,[1 1]*meanSampCurv,'k')

xlabel('Intensity Percentile')
ylabel([curvName ', 1/microns'])
%%Doesn't work because it has NaNs
%patch([curvX curvX(end:-1:1)],[curvCI(:,1)',curvCI(end:-1:1,2)'],'b')
legend(intName,intNameDN,'Mean');

plot(curvX,curvCIRandACN(:,1),'--b')
plot(curvX,curvCIRandACN(:,2),'--b')
plot(curvX,curvCIDNRandACN(:,1),'--r')
plot(curvX,curvCIDNRandACN(:,2),'--r')

xlim([min(curvX) max(curvX)])




figName = ['Artifactual corr and depth norm ' curvName ' vs  autocorr random plus noise nSig ' num2str(nSig) ' ac sig ' num2str(acSig) ' eSig ' num2str(eSig)];
if saveFigs,mfFigureExport(currFig,[outDir filesep figName]);end     

%% Example image from autocorrelated + noise random test image

currFig = figure;
tmp = testImRandACN;
%For visualization we set the voxels outside the cell to zero (this doesn't
%matter for the correlation because they aren't being sampled anyways)
tmp(~m) = 0;
imshow(tmp(:,:,iSlice),[])
figName = ['example plane autocorrelated random test image nSig ' num2str(nSig) ' ac sig ' num2str(acSig) ' eSig ' num2str(eSig)];

if saveFigs,mfFigureExport(currFig,[outDir filesep figName]);end    

%% depth varying image with spatial autocorrelation and noise

nRep = 10;
bpRandDepthACN = cell(nRep,1);
eSig = nSig/2;%Sigma of IID noise added on top
acSig = 2;
for iRep = 1:nRep
    
    tic
    % Depth - only variation image 
    testImDepthACN = bwdist(~m);
    nSig = 1;
    gSig = 4;
    amp = 10;
    testImDepthACN = amp*exp(-((testImDepthACN - gSig).^2) / gSig) + nSig*randn(size(testImDepthACN));
    testImDepthACN(m) = testImDepthACN(m) + amp/2;%Add an offset in the cell interior to avoid dividing by small numbers in the depth normalization
    %Spatially autocorrelated random image
    testImRandNoiseACN = randn(size(m));    
    testImRandNoiseACN = filterGauss3D(testImRandNoiseACN,acSig);
    %And set mean and std equal to the IID random image
    testImRandNoiseACN = testImRandNoiseACN - mean(testImRandNoiseACN(:));
    testImRandNoiseACN = testImRandNoiseACN ./ std(testImRandNoiseACN(:)) * nSig;
    testImRandNoiseACN = testImRandNoiseACN + amp/2;
    testImRandNoiseACN = testImRandNoiseACN + randn(size(m)) * eSig;
    testImDepthACN = testImDepthACN + testImRandNoiseACN;
    
    bpRandDepthACN{iRep} = analyze3DImageMaskedIntensities(testImDepthACN,m,sk,mp,2e3/MD.pixelSize_,false(size(m)));
            
    disp(iRep);
    toc
    
end

save([outDir filesep 'depth vary autocorr and noise test images and corr acSig ' num2str(acSig) ' eSig ' num2str(eSig) ' amp ' num2str(amp) '.mat'],'bpRandDepthACN','acSig','amp','nSig','eSig')
bpRandDepthACNAll = vertcat(bpRandDepthACN{:});

%% int vs. curv corr on depth varying and autocorrelated plus noise random test image

intField = 'intForCurvSampMean';
intName = 'Intensity';
intFieldDN = 'intForCurvSampMeanDepthNorm';
intNameDN = 'Depth Normalized Intensity';

%curvField = 'meanCurvSampMean';
curvField = 'MaxAbsPCCurvSampMean';
curvName = 'Max Absolute Principle Curvature';
curvConv = 1e3/MD.pixelSize_;



alpha = .05;
nBins = 100;
ptilesUse = linspace(0,100,nBins);

histBins = prctile(vertcat(bpRandDepthACNAll(:).(intField)),ptilesUse);
histBinsDN = prctile(vertcat(bpRandDepthACNAll(:).(intFieldDN)),ptilesUse);

meanCurvRandDepthACN= nan(nBins-1,1);
curvCIRandDepthACN = nan(nBins-1,2);
meanCurvDNRandDepthACN= nan(nBins-1,1);
curvCIDNRandDepthACN = nan(nBins-1,2);

meanCurvRandAllDepthACN= nan(nBins-1,nRep);
meanCurvDNRandAllDepthACN= nan(nBins-1,nRep);



for j = 1:(nBins-1)    
    
    for k = 1:nRep
        currSamp = bpRandDepthACNAll(k).(intField) > histBins(j) & bpRandDepthACNAll(k).(intField) <= histBins(j+1);    
        meanCurvRandAllDepthACN(j,k) = mean(bpRandDepthACNAll(k).(curvField)(currSamp) * curvConv);        

        currSampDN = bpRandDepthACNAll(k).(intFieldDN) > histBinsDN(j) & bpRandDepthACNAll(k).(intFieldDN) <= histBinsDN(j+1);    
        meanCurvDNRandAllDepthACN(j,k) = mean(bpRandDepthACNAll(k).(curvField)(currSampDN) * curvConv);        
    end
    
    meanCurvRandDepthACN(j) = mean(meanCurvRandAllDepthACN(j,:));
    bootSamp =  bootstrp(1e3,@mean, meanCurvRandAllDepthACN(j,:));
    curvCIRandDepthACN(j,:) = prctile(bootSamp,[alpha/2 100-alpha/2]);
    
    meanCurvDNRandDepthACN(j) = mean(meanCurvDNRandAllDepthACN(j,:));
    bootSamp =  bootstrp(1e3,@mean, meanCurvDNRandAllDepthACN(j,:));
    curvCIDNRandDepthACN(j,:) = prctile(bootSamp,[alpha/2 100-alpha/2]);
    
end

%% ---- Artifactual int vs. curv curves depth varying and autocorr plus noise --- %%

currFig = figure;
curvX = ptilesUse(1:end-1);
plot(curvX,meanCurvRandDepthACN)
hold on
plot(curvX,meanCurvDNRandDepthACN,'r')
meanSampCurv = nanmean(vertcat(bpRandDepthACNAll(:).(curvField)) * curvConv);
plot(xlim,[1 1]*meanSampCurv,'k')

xlabel('Intensity Percentile')
ylabel([curvName ', 1/microns'])
%%Doesn't work because it has NaNs
%patch([curvX curvX(end:-1:1)],[curvCI(:,1)',curvCI(end:-1:1,2)'],'b')
legend(intName,intNameDN,'Mean');

plot(curvX,curvCIRandDepthACN(:,1),'--b')
plot(curvX,curvCIRandDepthACN(:,2),'--b')
plot(curvX,curvCIDNRandDepthACN(:,1),'--r')
plot(curvX,curvCIDNRandDepthACN(:,2),'--r')

xlim([min(curvX) max(curvX)])




figName = ['Artifactual corr and depth norm ' curvName ' vs  depth var autocor rand noise nSig ' num2str(nSig) ' ac sig ' num2str(acSig) ' eSig ' num2str(eSig) ' amp ' num2str(amp)];
if saveFigs,mfFigureExport(currFig,[outDir filesep figName]);end     

%% Example image from depth varying autocorrelated + noise random test image

currFig = figure;
tmp = testImDepthACN;
%For visualization we set the voxels outside the cell to zero (this doesn't
%matter for the correlation because they aren't being sampled anyways)
tmp(~m) = 0;
imshow(tmp(:,:,iSlice),[])
figName = ['example plane depth var autocorr noise im nSig ' num2str(nSig) ' ac sig ' num2str(acSig) ' eSig ' num2str(eSig) ' amp ' num2str(amp)];

if saveFigs,mfFigureExport(currFig,[outDir filesep figName]);end    

%% ------ Intensity vs. Distance From Membrane ------ %%


outDir = 'W:\Hunter\orchestra_files_and_backup_merged\nih\Figures\Supplemental Figure Panels\Segmentation\Validation';

movPaths = 'W:\Hunter\orchestra_files_and_backup_merged\nih\4D gfpMIIB fix and stain\control\control10\movieData.mat';

MD = MovieData.load(movPaths,0);

iaDir = MD.processes_{MD.getProcessIndex('MaskedIntensity3DProcess',1,0)}.funParams_.OutputDirectory;

intAn = MD.processes_{MD.getProcessIndex('MaskedIntensity3DProcess',1,0)}.loadChannelOutput(1);

%%

chanInd = 1:2;%Just show the membrane channel so we don't have to deal with normalization and the apparent shifts it induces in these curves.
chanNames = {'MIIGFP','Membrane','Phalloidin'};
nChanShow = numel(chanInd);



for j = 1:nChanShow
    imNames = MD.channels_(chanInd(j)).getImageFileNames;
    imDir = MD.channels_(chanInd(j)).channelPath_;
    allIms(:,:,:,j) = stackRead([imDir filesep imNames{1}]);
    
    [backMean(j),backSTD(j)] = robustMean(double(reshape(allIms(:,:,:,j),1,[])),[],2);
    
end
    
    
    


%%
currFig = figure('Position',[ 549   -65   712   748]);
hold on
chanCols = [0 1 0 ; 
            0 0 1 ;
            1 0 0 ];
distVals = double(intAn.branchProfiles.wholeMaskDists);
nSpline = 1e3;
spDist = linspace(min(distVals),max(distVals),nSpline);        
distVals = distVals .* MD.pixelSize_ / 1e3;
distUnits = 'microns';

for j = 1:nChanShow
    currInt = intAn.branchProfiles.wholeMaskMeanVsDist(chanInd(j),:);
    currSTD = intAn.branchProfiles.wholeMaskSTDVsDist(chanInd(j),:);
    currInt = currInt ./ backMean(j);

    hold on
    plot(distVals,currInt,'color',chanCols(j,:),plotPars{:});
    
    xlabel(['Distance from Segmented Boundary, ' distUnits],labProps{:})
    ylabel('Fluroescence Intensity, Fold Over Background',labProps{:})
    
end


legend(chanNames(chanInd),labProps{:})

% for j = 1:nChanShow
%     
%     currInt = intAn.branchProfiles.wholeMaskMeanVsDist(chanInd(j),:);
%     currSTD = intAn.branchProfiles.wholeMaskSTDVsDist(chanInd(j),:);
%     plot(distVals,currInt + currSTD,'--','color',chanCols(j,:));
%     plot(distVals,currInt - currSTD,'--','color',chanCols(j,:));    
% end

xlim([min(distVals) max(distVals)])

currAx = get(currFig,'CurrentAxes');
set(currAx,axProps{:})
plot([0 0],ylim,'k',plotPars{:})
figName = ['intensity vs distance from segmented boundary example'];
if saveFigs,mfFigureExport(currFig,[outDir filesep figName]);end    





