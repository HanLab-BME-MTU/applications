function [] = makeMovieProtrusionOverlay(MO)
% This is a wrapper function for makeProtrusionOverlay
if isa(MO,'MovieData')
    MO = MovieData.load(MO.getFullPath);
    MDAll{1} = MO;
    nMovies = numel(MDAll);
elseif isa(MO,'MovieList')
    ML = MovieList.load(MO.getFullPath);
    MDAll = ML.movies_;
    nMovies = numel(MDAll);
end
%% ROI for cell to reduce memory burden
for ii=1:nMovies
    MD = MDAll{ii};

    % Get the windowing package
    windPack = MD.getPackage(MD.getPackageIndex('WindowingPackage'));
    % Protrusion sampling process
    protProc = windPack.processes_{1};
    % Load smoothedEdge
    smoothedEdge = protProc.loadChannelOutput('output','smoothedEdge');
    % Load the mask refinement process
    maskProc = MD.getProcess(MD.getProcessIndex('MaskRefinementProcess'));
    % Load first img
    iChan = find(maskProc.checkChannelOutput);
    img = MD.channels_(iChan).loadImage(1);
    minIm = double(quantile(img(:),0.04));
    maxIm = double(quantile(img(:),0.998));

    figPath= [MD.getPath filesep 'EdgeOverlay']; mkdir(figPath)
    
    h=makeProtrusionOverlay(smoothedEdge,img,minIm,maxIm,MD.timeInterval_);
    
    print(h, '-depsc2', strcat(figPath,'/coloredEdgeOverlay.eps'));
    print(h, '-dtiff', strcat(figPath,'/coloredEdgeOverlay','.tif'));
    savefig(h,strcat(figPath,'/coloredEdgeOverlay.fig'),'compact')
end

