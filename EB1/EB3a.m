function EB3a


for i = 1:1

    s = 3;
    strg=sprintf('%%.%dd',s);
    indxStr=sprintf(strg,i);

    %     I = imread(['Z:\AlexData\EB1_Katsu\20041113EB1-GFP',indxStr,'.tif']);
    %     I = imread(['Z:\AlexData\EB1_Claudio\controlmovie5\EB1-GFP_control05_T',indxStr,'.tif']);
    %     I = imread(['Z:\AlexData\EB1_Lisa\extract\extract_tiffs\xjul03_r11\xjul03_r11tiff',indxStr,'.tif']);

    %        I = imread(['Z:\AlexData\EB1_Katsu_2s_1\EB1-GFP-5',indxStr,'.tif']);
    %           I = imread(['Z:\AlexData\EB1_Katsu_2s_2\EB1-GFP-6',indxStr,'.tif']);
    %-----CLAUDIO--------------
%     I = imread(['C:\amatov\data\Monastrol\images\mono3',indxStr,'.tif']);
        I = imread(['C:\amatov\data\070622_2\images\mono1',indxStr,'.tif']);
    % I = imread(['Z:\AlexData\Ben\EB1tracking\Op18Spindles\102406_Op18\Op18_movie_2\images\Op18_movie2_',indxStr,'.tif']);
%         I = imread(['X:\AlexData\Yukako\060707_2\2_',indxStr,'.tif']);
%     I = imread(['X:\AlexData\Jay\EB1_60x2x2binning_25uMSTC_4a\images\EB1_60x_2x2binning_25uMSTC_4a',indxStr,'.tif']);
    
%     I = imread(['P:\forAlex\images300\RNAi_RACCA_3_',indxStr,'.tif']);
  
% I = imread(['X:\AlexData11\786Opar\786Opar_NaCl01_R3D\images\786Opar_NaCl01_T',indxStr,'.tif']);
    %------KATSU----------------
    % I = imread(['Z:\AlexData\Katsu\cell7\images\EB1-GFP',indxStr,'.tif']);

    I=double(I);
    aux = Gauss2D(I,1);%1 
    I2 = Gauss2D(I,4); %4 (Yukako 10)
    I3 = aux - I2;
    [cutoffInd, cutoffV] = cutFirstHistMode(I3,0);

    coef = 1;% coef = 4 Katsu; coef = 1 Claudio; coef = 1 Lisa_xju103_r11; 
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

    %-------------------------------------
    %     aaux = 5;
    %     If=Gauss2D(I,1);
    %     figure, imshow(If(1+aaux:end-aaux,1+aaux:end-aaux),[]);%I4
    %     hold on
    %     for i = 1:length(feats.ori)
    %         h = quiver(feats.pos(i,1)-aaux,feats.pos(i,2)-aaux,-cos(feats.ori(i)*pi/180),sin(feats.ori(i)*pi/180),3,'r');
    %         set(h,'LineWidth',2)
    %     end
    %------------------------------------

    Cm = nanmean(Crop,3); % MEAN/REPRESENTATIVE EB1 CROP
    Crop(isnan(Crop))=0;% border effect - some NaN
    Cm1 = bwlabel(Cm);
    statsC = regionprops(Cm1,'all');

    %     figure,imshow(Cm,[])
    %
    %     uiviewpanel,imshow(Cm,[])


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

    %     featNames = fieldnames(feats);
    %     for field = 1:length(featNames)
    %         feats.(featNames{field}) = feats.(featNames{field})(goodFeats,:);
    %     end

    %     CC = Crop(:,:,goodFeats);
    %     CCm = nanmean(CC,3); %MEAN EB1 BASED on GOOD/FITTED EB1 SPOTS
    %     CC(isnan(CC))=0;% border effect - some NaN
    %     CCm1 = bwlabel(CCm);
    %     statsCC = regionprops(CCm1,'all');

    %     figure,imshow(CCm,[])
    %
    %     figure,histogram(ssqC)
    %     figure,histogram(ssqC(goodFeats))
    %     [cutoffIndex, cutoffValue] = cutFirstHistMode(ssqC(goodFeats),0);
    %     [cutoffIndex, cutoffValue] = cutFirstHistMode(ssqC,0);
    %
    % %     veryGoodFeats = find(ssqC(goodFeats)<cutoffValue);
    %     veryGoodFeats = find(ssqC<cutoffValue);



    featNames = fieldnames(feats);
    for field = 1:length(featNames)
        feats.(featNames{field}) = feats.(featNames{field})(goodFeats,:);
    end



    %-------------------------------------
    aaux = 5;
    If=Gauss2D(I,1);
    figure, imshow(If(1+aaux:end-aaux,1+aaux:end-aaux),[]);%I4
    hold on
    for i = 1:length(feats.ori)
        h = quiver(feats.pos(i,1)-aaux,feats.pos(i,2)-aaux,-cos(feats.ori(i)*pi/180),sin(feats.ori(i)*pi/180),3,'r');
        set(h,'LineWidth',2)
    end
    %------------------------------------
    %        save(['Z:\AlexData\Torsten\control_1\cands\feats',indxStr],'feats')
%     save(['C:\amatov\data\Monastrol\cands\feats',indxStr],'feats')
%     
%     clear all
%     close all
end
% p
% aaux = 5;
% If=Gauss2D(I,1);
% figure, imshow(If(1+aaux:end-aaux,1+aaux:end-aaux),[]);%I4
% hold on
% for i = 1:length(feats.ori)
%     h = quiver(feats.pos(i,1)-aaux,feats.pos(i,2)-aaux,-cos(feats.ori(i)*pi/180),sin(feats.ori(i)*pi/180),3,'r');
%     set(h,'LineWidth',2)
% end

% stats
%
%
%
% figure,imshow(I,[])
% hold on
% plot(xGrid(:),yGrid(:),'g.')
%
% figure,imshow(Crop1,[])


