function EB3a(debug,coef,sigma)

% input debug -> if 0 - save 125 cands files to disk, else if 1 - run for
% the 1st file and plot detection figure
% 
% coef = 1; sigma = 4;
%
% run:
% EB3a(1,1,4) to plot figure and look at results
% EB3a(0,1,4) to run thru whole movie and save detection

[fileName,dirName] = uigetfile('*.tif','Choose a .tif file');
I = imread([dirName,filesep,fileName]);
if debug == 0
    le = 125
elseif debug == 1
    le = 1
end
for i = 1:le

    s = 3;
    strg=sprintf('%%.%dd',s);
    indxStr=sprintf(strg,i);

    I = imread([dirName,fileName(1:end-7),indxStr,'.tif']);
    I=double(I);
    aux = Gauss2D(I,1);%1 
    I2 = Gauss2D(I,sigma); %4 (Yukako 10)
    I3 = aux - I2;
    [cutoffInd, cutoffV] = cutFirstHistMode(I3,0);

    % coef = 4 Katsu; coef = 1 Claudio; coef = 1 Lisa_xju103_r11; 
    I4 = I3>cutoffV*coef; % REMOVE THE NOISE FEATURES %no 3

    X = bwlabel(I4);
    warningState = warning;
    warning off all
    stats = regionprops(X,'all'); % Warning: Out of range value converted to intmin('uint8') or intmax('uint8').
    warning(warningState)

    % Initialize 'feats' structure
    feats=struct(...
        'pos',[0 0],...                  % Centroid - [y x]
        'ecc',0,...                      % Eccentricity
        'ori',0);   % Orientation

    for j = 1:length(stats)
        feats.pos(j,1) = stats(j).Centroid(1);
        feats.pos(j,2) = stats(j).Centroid(2);
        feats.ecc(j,1) = stats(j).Eccentricity;
        feats.ori(j,1) = stats(j).Orientation;
        feats.len(j,1) = stats(j).MajorAxisLength;

        e1 = [-cos(stats(j).Orientation*pi/180) sin(stats(j).Orientation*pi/180) 0];
        e2 = [sin(stats(j).Orientation*pi/180) cos(stats(j).Orientation*pi/180) 0];
        e3 = [0 0 1];
        Ori = [stats(j).Centroid  0];
        v1 = [-10 10];
        v2 = [-5 5];
        v3 = [0 0];
        [xGrid,yGrid]=arbitraryGrid(e1,e2,e3,Ori,v1,v2,v3);

        Crop(:,:,j) = interp2(I,xGrid,yGrid);
        %         Crop(:,:,j) = interp2(I,xGrid,yGrid,'*linear');

        e1 = [];e2 = [];e3 = []; Ori = []; v1 = []; v2 = []; xGrid = []; yGrid = [];
    end

    Cm = nanmean(Crop,3); % MEAN/REPRESENTATIVE EB1 CROP
    Crop(isnan(Crop))=0;% border effect - some NaN
    Cm1 = bwlabel(Cm);
    statsC = regionprops(Cm1,'all');

    sC = size(Crop);
    Cm3d = repmat(Cm,[1,1,size(Crop,3)]);
    dC = Crop - Cm3d;
    sqC = dC.^2;
    ssqC = squeeze(sum(sum(sqC,1),2)); %LIST OF DIFFERENCES AFTER SUBTRACTION

    B = Cm(:);
    A = ones(length(B),2);

    for m = 1:size(Crop,3)
        CR = Crop(:,:,m);
        A(:,2) = CR(:);
        goodRows = find(A(:,2) ~= 0 & isfinite(B));
        XX = lscov(A(goodRows,:),B(goodRows));
        RES = B(goodRows) - A(goodRows,:)*XX;
        OUT(m,:) = [mean(RES(:).^2),XX'];
    end

    [Ind,V]=cutFirstHistMode(OUT(:,1),0);% switch to 1 to see HIST

    goodFeats = find(OUT(:,1)<V); % SPOTS WHICH FIT WELL WITH THE MEAN EB1 SPOT

    featNames = fieldnames(feats);
    for field = 1:length(featNames)
        feats.(featNames{field}) = feats.(featNames{field})(goodFeats,:);
    end

    if debug == 1
        aaux = 5;
        If=Gauss2D(I,1);
        figure, imshow(If(1+aaux:end-aaux,1+aaux:end-aaux),[]);%I4
        hold on
        for i = 1:length(feats.ori)
            h = quiver(feats.pos(i,1)-aaux,feats.pos(i,2)-aaux,-cos(feats.ori(i)*pi/180),sin(feats.ori(i)*pi/180),3,'r');
            set(h,'LineWidth',2)
        end
    elseif debug == 0
        save([dirName,'cands\feats',indxStr],'feats')
        clear all
        close all
    end
end


