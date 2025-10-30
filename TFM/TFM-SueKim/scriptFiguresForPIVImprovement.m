% Sue's code for saving all publishable figures

%% Displacement field for our original PTV
% clc
% clear all

displOriginalPTVstruct = load('/storage/network/TFM_Development/TFM2D/PIVimprovement/simulation/TFMPackage/displacementField/displField.mat');
displOriginalPTV = displOriginalPTVstruct.displField;

% quiver plot - original PTV
figure(1), quiver(displOriginalPTV(1).pos(:,1),displOriginalPTV(1).pos(:,2),...
    displOriginalPTV(1).vec(:,1),displOriginalPTV(1).vec(:,2),0,'Color','g')
hold on

% missing spots
j=1;
unTrackedBeads=isnan(displOriginalPTV(j).vec(:,1));
% sumUntracked = sum(unTrackedBeads);
currentBeads = displOriginalPTV(j).pos(unTrackedBeads,:);
plot(currentBeads(:,1),currentBeads(:,2),'ro')
hold on

% Filtered vectors
displOriginalPTVFilteredstruct = load('/storage/network/TFM_Development/TFM2D/PIVimprovement/simulation/TFMPackage/correctedDisplacementField/displField.mat');
displOriginalPTVFiltered = displOriginalPTVFilteredstruct.displField;

trackedBeadsBeforeFiltering = ~unTrackedBeads;
trackedBeadsAfterFiltering=~isnan(displOriginalPTVFiltered(j).vec(:,1));
filteredBeads = trackedBeadsBeforeFiltering & ~trackedBeadsAfterFiltering;

quiver(displOriginalPTVFiltered(j).pos(trackedBeadsAfterFiltering,1),displOriginalPTVFiltered(j).pos(trackedBeadsAfterFiltering,2),...
    displOriginalPTVFiltered(j).vec(trackedBeadsAfterFiltering,1),displOriginalPTVFiltered(j).vec(trackedBeadsAfterFiltering,2),0,'Color','b')

%% Displacement map for our original PTV
dataPath='/storage/network/TFM_Development/TFM2D/PIVimprovement/simulation/Figures/displMapOriginalPTVStep2';
generateHeatmapFromField(displOriginalPTV,dataPath,0,30); % step 2
dataPath='/storage/network/TFM_Development/TFM2D/PIVimprovement/simulation/Figures/displMapOriginalPTVStep3';
generateHeatmapFromField(displOriginalPTVFiltered,dataPath,0,30); % step 
%% Displacement field produced from ImageJ
% Import .txt file of ImageJ 
x = fileread('ImageJ_PIV_Stack_final.txt');
str = convertCharsToStrings(x);
str_sp = split(str);

% Extract X and Y Values
X = str_sp(1:16:end-1);
Y = str_sp(2:16:end);

ux_ImageJ = str_sp(3:16:end);
uy_ImageJ = str_sp(4:16:end);

X_convert = {X};
S = sprintf('%s ', X_convert{:});
X = sscanf(S, '%f');

Y_convert = {Y};
S = sprintf('%s ', Y_convert{:});
Y = sscanf(S, '%f');

ux_convert = {ux_ImageJ};
S = sprintf('%s ', ux_convert{:});
ux_ImageJ = sscanf(S, '%f');

uy_convert = {uy_ImageJ};
S = sprintf('%s ', uy_convert{:});
uy_ImageJ = sscanf(S, '%f');

% Displacement Field ImageJ Result
positionImageJ = [X,Y];
vectorImageJ = [ux_ImageJ,uy_ImageJ];
% Save
% save('dispFieldImageJ.mat','dispFieldImageJ')
save('dispPositionImageJ.mat','positionImageJ')
save('dispVectorImageJ.mat','vectorImageJ')

%% Displacement filed for ImageJ PIV Results
% displOriginalImageJstruct = load('/storage/network/TFM_Development/TFM2D/PIVimprovement/simulation/TFMPackage/displacementFieldImageJ/dispFieldImageJ.mat');
dispPositionImageJstruct = load('dispPositionImageJ.mat');
dispVectorImageJstruct = load('dispVectorImageJ.mat');

%% quiver plot - ImageJ result
figure(2), quiver(dispPositionImageJstruct.positionImageJ(:,1),dispPositionImageJstruct.positionImageJ(:,2),...
    dispVectorImageJstruct.vectorImageJ(:,1),dispVectorImageJstruct.vectorImageJ(:,2),0,'Color','b')
hold on
%% color map - ImageJ
% generateHeatmapFromGridData(positionImageJ(:,1),positionImageJ(:,2),vectorImageJ(:,1),vectorImageJ(:,2),0,0,6000,false,430,430);
dataPath='/storage/network/TFM_Development/TFM2D/PIVimprovement/simulation';
generateHeatmapFromGridData(dispPositionImageJstruct.positionImageJ(:,1),dispPositionImageJstruct.positionImageJ(:,2),dispVectorImageJstruct.vectorImageJ(:,1),dispVectorImageJstruct.vectorImageJ(:,2),[dataPath '/Figures'],0,0,30,false,500,500);
%% Compare between ground-truth displacement vectors and ImageJ output
% load ground-truth displacement vector field
uxStruct = load('ux_sim.mat');
ux = uxStruct.ux;
uyStruct = load('uy_sim.mat');
uy = uyStruct.uy;

% for loop to find similar vector values with tolerance
% the size of ground truth displacement is 512x512
% the size of PIV (ImageJ) is 900x2 for pos and vec (dispPositionImageJ and dispVectorImageJ

% nElements = length(dispVectorImageJstruct.vectorImageJ);   %900
% for i = 1:nElements
%     matches = abs(displOriginalPTV.vec - dispVectorImageJstruct.vectorImageJ(i)) < tol;
% end

%% 1. Find out positions that miss the vectors
emptyLocations = isnan(dispVectorImageJstruct.vectorImageJ(:,1)) | isempty(dispVectorImageJstruct.vectorImageJ(:,1));
% sum(empytyLocations) = 0
%% 2. Find out vectors that are out of tolerance.
umag = (ux.^2+uy.^2).^0.5;
% tol = std(umag);     %tolerance threshold 70-80%
% nElements = length(dispVectorImageJstruct.vectorImageJ);   %900
% uDownSampled = zeros(size(dispVectorImageJstruct.vectorImageJ));
% matchingVector = true(size(dispVectorImageJstruct.vectorImageJ(:,1)));
% for ii = 1:nElements
%     % Find the location at which dispVectorImageJstruct.vectorImageJ(i)belong to
%     curPos = dispPositionImageJstruct.positionImageJ(ii,:);
%     curU = [ux(curPos(2)+1,curPos(1)+1),uy(curPos(2)+1,curPos(1)+1)];  % assign ground-truth disp. vector onto positionImageJ
%     curDiffVec = curU - dispVectorImageJstruct.vectorImageJ(ii);    %ground-truth vector - ImageJvector
%     matchingVector(ii) = (curDiffVec(1)^2+curDiffVec(2)^2)^0.5 < tol;%*(curU(1)^2+curU(2)^2)^0.5;
% %     matches_large = (curDiffVec(1)^2+curDiffVec(2)^2)^0.5 > tol*(curU(1)^2+curU(2)^2)^0.5;
%     uDownSampled(ii,:)=curU;
% end
% unMatching = ~matchingVector;

tol = 0.50;     %tolerance threshold 70-80%
nElements = length(dispVectorImageJstruct.vectorImageJ);   %900
uDownSampled = zeros(size(dispVectorImageJstruct.vectorImageJ));
matchingVector = true(size(dispVectorImageJstruct.vectorImageJ(:,1)));
for ii = 1:nElements
    % Find the location at which dispVectorImageJstruct.vectorImageJ(i)belong to
    curPos = dispPositionImageJstruct.positionImageJ(ii,:);
    curU = [ux(curPos(2),curPos(1)),uy(curPos(2),curPos(1))];  % assign ground-truth disp. vector onto positionImageJ
    curDiffVec = curU - dispVectorImageJstruct.vectorImageJ(ii);    %ground-truth vector - ImageJvector
    matchingVector(ii) = (curDiffVec(1)^2+curDiffVec(2)^2)^0.5 < tol*(curU(1)^2+curU(2)^2)^0.5;  %magnitude
    %     matches_large = (curDiffVec(1)^2+curDiffVec(2)^2)^0.5 > tol*(curU(1)^2+curU(2)^2)^0.5;
    uDownSampled(ii,:)=curU;
end
unMatching = ~matchingVector;
quiver(dispPositionImageJstruct.positionImageJ(unMatching,1),dispPositionImageJstruct.positionImageJ(unMatching,2),...
dispVectorImageJstruct.vectorImageJ(unMatching,1),dispVectorImageJstruct.vectorImageJ(unMatching,2),0,'Color','g')
hold on
%% 3. Plot the missing positions as a circle
plot(dispPositionImageJstruct.positionImageJ(emptyLocations,1),dispPositionImageJstruct.positionImageJ(emptyLocations,2),'ro')
%% 4. Plot the wrong vectors with a red color
quiver(dispPositionImageJstruct.positionImageJ(unMatching,1),dispPositionImageJstruct.positionImageJ(unMatching,2),...
    dispVectorImageJstruct.vectorImageJ(unMatching,1),dispVectorImageJstruct.vectorImageJ(unMatching,2),0,'Color','r')
plot(dispPositionImageJstruct.positionImageJ(unMatching,1),dispPositionImageJstruct.positionImageJ(unMatching,2),'ro')

%% quiver plot of filtered vectors
quiver(dispPositionImageJstruct.positionImageJ(:,1),dispPositionImageJstruct.positionImageJ(:,2),...
    uDownSampled(:,1), uDownSampled(:,2),0,'Color','r')
hold on

% missing spots
% %Find elements in one array not in another
% nElements = length(dispVectorImageJstruct.vectorImageJ);   %900
% for i = 1:nElements
%     filteredVectors = abs(uDownSampled(:,1) - dispVectorImageJstruct.vectorImageJ(i,1)) > tol;    
% end
% 
% % currentBeadsJ = dispPositionImageJstruct(j).positionImageJ(filteredVectors,:);
% % plot(currentBeadsJ(:,1),currentBeadsJ(:,2),'ro')
% 
% % sum(filteredVectors) = 846
% j=1;
% filtered_wrong = dispVectorImageJstruct(j).vectorImageJ(~filteredVectors,:);
% filtered_correct = dispVectorImageJstruct(j).vectorImageJ(filteredVectors,:);
% quiver(dispPositionImageJstruct(j).positionImageJ(filtered_correct,1),dispPositionImageJstruct(j).positionImageJ(filtered_correct,2),...
%     dispVectorImageJstruct(j).vectorImageJ(filtered_correct,1),dispVectorImageJstruct(j).vectorImageJ(filtered_correct,2),0,'Color',[0.5 0.5 0.5])
% plot(dispPositionImageJstruct.positionImageJ(unMatching,1),dispPositionImageJstruct.positionImageJ(unMatching,2),'ro')

% displCorrectedPTVFilteredstruct = load('/storage/network/TFM_Development/TFM2D/PIVimprovement/simulation/TFMPackage/correctedDisplacementField/displField.mat');
% displCorrectedPTVFiltered = displCorrectedPTVFilteredstruct.displField;
% 
% unTrackedBeads_J=isnan(displCorrectedPTVFiltered(j).vec(:,1));
% trackedBeadsBeforeFiltering_J = ~unTrackedBeads_J;
% trackedBeadsAfterFiltering_J = ~isnan(displCorrectedPTVFiltered(j).vec(:,1));
% filteredBeads_J = trackedBeadsBeforeFiltering_J & ~trackedBeadsAfterFiltering_J;
% 
% % trackedBeadsAfterFiltering_J=displCorrectedPTVFiltered(j).vec(:,1);
% 
% quiver(displCorrectedPTVFiltered(j).pos(trackedBeadsAfterFiltering_J,1),displCorrectedPTVFiltered(j).pos(trackedBeadsAfterFiltering_J,2),...
%     displCorrectedPTVFiltered(j).vec(trackedBeadsAfterFiltering_J,1),displCorrectedPTVFiltered(j).vec(trackedBeadsAfterFiltering_J,2),0,'Color',[0.5 0.5 0.5])

%% 


matches = (uDownSampled(1)^2+uDownSampled(2)^2)^0.5 < tol*(uDownSampled(1)^2+uDownSampled(2)^2)^0.5;  %magnitude
j=1;
%matches = untrackedbeads
currentBeads = disp

% trackedBeads_J = positionImageJ(mathces);
trackedBeads_VecJ = uDownSampled(matches,:);
unTrackedBeads_J = dispPositionImageJstruct(j).positionImageJ(matches)

trackedBeads_VecJ = positionImageJ(matches,:);

plot(trackedBeads_VecJ(:,1),trackedBeads_VecJ(:,2),'ro')
%%
c = setdiff(dispVectorImageJstruct.vectorImageJ,uDownSampled);

[LIA] = ismembertol(dispVectorImageJstruct.vectorImageJ(:,1),uDownSampled(:,1),0.2,'OutputAllIndices',true);

j=1;
currentBeadsPos_J = dispPositionImageJstruct(j).positionImageJ(LIA,:); 
currentBeadsVec_J = dispVectorImageJstruct(j).vectorImageJ(LIA,:);

% j =1;
% untrackedBeadsPos_J = dispPositionImageJstruct(j).positionImageJ(~matches,:);

c = intersect(dispVectorImageJstruct.vectorImageJ,uDownSampled,'rows');

c1 = setdiff(dispVectorImageJstruct.vectorImageJ(:,1),uDownSampled(:,1));
c2 = setdiff(vectorImageJ(:,2),uDownSampled(:,2));
unTrackedBeads_J = vectorImageJ(~uDownSampled);



%sumMatches = sum(matches) % 145
%vectorImageJ 900x2 
%uDownSampled 900x2


unTrackedBeads_J=dispVectorImageJstruct.vectorImageJ(~uDownSampled,:);


currentBeadsPos_J = displOriginalPTV(j).pos(matches,:);  
currentBeadsVec_J = displOriginalPTV(j).vec(matches,:); 
figure, quiver(currentBeadsPos_J(:,1),currentBeadsPos_J(:,2),...
    currentBeadsVec_J(:,1),currentBeadsVec_J(:,2),0,'Color',[0.5 0.5 0.5])
hold on
untrackedBeadsPos_J = displOriginalPTV(j).pos(~matches,:);
untrackedBeadsVec_J = displOriginalPTV(j).vec(~matches,:);

quiver(untrackedBeadsPos_J(:,1),untrackedBeadsPos_J(:,2),...
    untrackedBeadsVec_J(:,1),untrackedBeadsVec_J(:,2),0,'Color','r')


%% mpiv
%load displacement field
iuStruct = load('iu.mat'); % displacement vector calculated by mpiv
ivStruct = load('iv.mat');  % displacement vector calculated by mpiv
iu = iuStruct.iu;
iv = ivStruct.iv;

%load position
xiStruct = load('xi.mat'); % 1x1 struct
yiStruct = load('yi.mat'); % 1x1 struct
xi = xiStruct.xi;
yi = yiStruct.yi;

%make them 30x30
[xgrid,ygrid] = meshgrid(xi,yi);

%quiver plot
figure, quiver(xgrid,ygrid,iu,iv,0,'Color','g')
hold on
%% mpiv - load filtered displacement field
iu_ft_Struct = load('iu_ft.mat');
iv_ft_Struct = load('iv_ft.mat');
iu_ft = iu_ft_Struct.iu_ft;
iv_ft = iv_ft_Struct.iv_ft;

quiver(xgrid,ygrid,iu_ft,iv_ft,0,'Color','b')
hold on
%% mpiv - load interpolated displacement field
load('missingPos_mpiv.mat')
iu_ip_Struct = load('iu_ip.mat');
iv_ip_Struct = load('iv_ip.mat');
iu_ip = iu_ip_Struct.iu_ip;
iv_ip = iv_ip_Struct.iv_ip;

plot(missingPos_mpiv(:,1),missingPos_mpiv(:,2),'ro')
% % x = load('missingPosVec_mpiv.mat');
% % missingVec = x.missingVec_mpiv;
% % currentB = x.missingPos_mpiv(missingVec_mpiv,:);
% % figure, plot(currentB(:,1),currentB(:,2),'ro')

% names = fieldnames(x);
% xx = x.(names{:});

% x = uiimport('test_signal.mat');
% names =fieldnames(x);
% xx = x.(names{:});
% total_input_signal = [zeros(1,size(xx,2));xx];
% 
% 
% figure, quiver(missingPos_mpiv(:,1),missingPos_mpiv(:,2),missingVec_mpiv(:,1),missingVec_mpiv(:,2),0,'ro')

% mpivVec = [iu_ip iv_ip];
% mpivPos = [xgrid ygrid];
% j = 1;
% missingBeads=isnan(mpivVec(:,1));
% current = xgrid(missingBeads,:);
% plot(mpivPos(:,1),mpivPos(:,2),'ro');

% quiver(missingPos_mpiv(:,1),missingPos_mpiv(:,2),missingVec_mpiv(:,1),missingVec_mpiv(:,2),0,'ro')
% c = setdiff(iu_ip,iu_ft);
% d = setdiff(iv_ip, iv_ft);
% 
% missingVec = [c d];
% 
% 
% plot(c,d,'ro')

% quiver(xgrid,ygrid,iu_ip,iv_ip,0,'Color','r')
% hold on
% %% mpiv - load filtered displacement vector
% u_mpiv_filtered = load('iu_ft.mat'); %filtered displacement vector by mpiv
% v_mpiv_filtered = load('iv_ft.mat'); %filtered displacement vector by mpiv
% 
% iu_ft = u_mpiv_filtered.iu_ft;
% iv_ft = v_mpiv_filtered.iv_ft;
% 
% quiver(xgrid,ygrid,iu_ft,iv_ft,0,'Color','m')


%% mpiv - get missing coordinates
% figure, imshow(ones(500))
% hold on
% posNA = [];
% n = 0;
% % Loop, picking up the points.
% disp('Left mouse button picks points.')
% disp('Right mouse button picks last point.')
% 
% but = 1;
% while but == 1
%     [xi,yi,but] = ginput(1);
%     plot(xi,yi,'ro')
%     n = n+1;
%     text(xi,yi-8,num2str(n));
%     posNA(n,:) = [xi yi];
% end
% hold off
% %% posFA- get the coordinates
% figure, imshow(ones(500))
% 
% hold on
% posFA = [];
% n = 0;
% % Loop, picking up the points.
% disp('Left mouse button picks points.')
% disp('Right mouse button picks last point.')
% 
% but = 1;
% while but == 1
%     [xi,yi,but] = ginput(1);
%     plot(xi,yi,'ro')
%     n = n+1;
%     text(xi,yi-8,num2str(n));
%     posFA(n,:) = [xi yi];
% end
% hold off
%% Figure 1(d) PIV quiver plot
%Run TFM before run this code. PIV Suite
load('/storage/network/TFM_Development/TFM2D/PIVimprovement/simulation/TFMPackage/displacementField/displField.mat')
figure, quiver(displField.pos(:,1),displField.pos(:,2),...
    displField.vec(:,1),displField.vec(:,2),0,'Color','b')
hold on
%% colormap - mpiv
dataPath='/storage/network/TFM_Development/TFM2D/PIVimprovement/simulation';
generateHeatmapFromGridData(xgrid,ygrid,iu_ft,iv_ft,[dataPath '/Figures'],0,0,30,false,430,430);
hold on
generateHeatmapFromGridData(xgrid,ygrid,iufilled,ivfilled,[dataPath '/Figures'],0,0,30,false,500,500);
%% PIV - load corrected displacement field
load('/storage/network/TFM_Development/TFM2D/PIVimprovement/simulation/TFMPackage/correctedDisplacementField/displField.mat')
figure, quiver(displField.pos(:,1),displField.pos(:,2),...
    displField.vec(:,1),displField.vec(:,2),0,'Color','b')
%% PIV - load calculated displacement field
load('/storage/network/TFM_Development/TFM2D/PIVimprovement/simulation/TFMPackage/displacementField/displField.mat')
figure, quiver(displField.pos(:,1),displField.pos(:,2),...
    displField.vec(:,1),displField.vec(:,2),0,'Color','b')