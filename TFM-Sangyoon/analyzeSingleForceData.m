function [meanDispError,dispDetec,meanForceError,peakForceRatio,forceDetec] = analyzeSingleForceData(f,d,minCorLength,dataPath)
%% single force experiment
% input parameters to be replaced with function inputs
% f=2000; %Pa
% d=10;
% minCorLength = 21;
imgPath=[dataPath filesep 'Beads'];
refPath=[dataPath filesep 'Reference'];
orgPath=[dataPath filesep 'Original'];
analysisFolder = dataPath;
if ~exist(refPath,'dir') || ~exist(orgPath,'dir')
    mkdir(imgPath);
    mkdir(refPath);
    mkdir(orgPath);
end
%% reference image (300x200)
xmax=200;
ymax=300;

nPoints = 5000;
% bead_r = 40; % nm
% pixSize = 72e-9; % nm/pix 90x
%% Now displacement field from given force
E=8000;  %Young's modulus, unit: Pa
meshPtsFwdSol=2^10;
forceType = 'groupForce';

gridSpacing = 1;
xmin = gridSpacing;
ymin = gridSpacing;

%% loading original displacement field and force field
load([orgPath filesep 'data.mat'],'ux','uy','x_mat_u','y_mat_u','bead_x','bead_ux','bead_y','bead_uy','force_x','force_y');

%% Now force reconstruction via movieData (non-GUI mode)
% Retrieve current location
MD=MovieData.load(fullfile(analysisFolder,'movieData.mat'));

%% Create TFM package and retrieve package index
iPack=  MD.getPackageIndex('TFMPackage');

%% Postprocessing - saving and analyzing force field
% Loading displacement field and force field
% Load the displField
disp('Calculating displacement errors and force errors...')
displField=MD.getPackage(iPack).getProcess(3).loadChannelOutput;
% finding displacement at bead location
org_ux = zeros(size(displField(1).pos(:,1)));
org_uy = zeros(size(displField(1).pos(:,1)));
nmPoints = length(displField(1).pos(:,1));
for k=1:nmPoints
    [~,indcol_closest_x] = min(abs(x_mat_u(1,:)-displField(1).pos(k,1)),[],2);
    [~,indrow_closest_y] = min(abs(y_mat_u(:,1)-displField(1).pos(k,2)),[],1);
    row_bottom = max(1,indrow_closest_y-3);
    row_top = min(size(x_mat_u,1),indrow_closest_y+3);
    col_bottom = max(1,indcol_closest_x-3);
    col_top = min(size(y_mat_u,2),indcol_closest_x+3);
    loc_xmat = x_mat_u(row_bottom:row_top,col_bottom:col_top);
    loc_ymat = y_mat_u(row_bottom:row_top,col_bottom:col_top);
    loc_ux = ux(row_bottom:row_top,col_bottom:col_top);
    loc_uy = uy(row_bottom:row_top,col_bottom:col_top);
    org_ux(k) = interp2(loc_xmat,loc_ymat,loc_ux,displField(1).pos(k,1),displField(1).pos(k,2));
    if isnan(org_ux(k))
        org_ux(k) = ux(indrow_closest_y,indcol_closest_x);
    end
    org_uy(k) = interp2(loc_xmat,loc_ymat,loc_uy,displField(1).pos(k,1),displField(1).pos(k,2));
    if isnan(org_uy(k))
        org_uy(k) = uy(indrow_closest_y,indcol_closest_x);
    end
end
%% errors in displacementfield
displField.vec(isnan(displField.vec(:,1)),:) = 0;
meanDispError= nansum(((org_ux-displField(1).vec(:,1)).^2+(org_uy-displField(1).vec(:,2)).^2).^.5)/(length(displField(1).vec(:,2))); %normalized by the number of beads% detectability (u at force application / u at background)
maskForce = ((x_mat_u-100).^2+(y_mat_u-150).^2)<=d/2*2;
dispDetecIdx = maskVectors(displField(1).pos(:,1),displField(1).pos(:,2),maskForce);
if isempty(dispDetecIdx)
    dispDetec = 0;
else
    displFieldForce = displField(1).vec(dispDetecIdx,:);
    displFieldMag = (displFieldForce(:,1).^2+displFieldForce(:,2).^2).^0.5;
    backgroundIdx = maskVectors(displField(1).pos(:,1),displField(1).pos(:,2),~bwmorph(maskForce,'dilate',floor(d/2)));
    displFieldBgd = displField(1).vec(backgroundIdx,:);
    displFieldBgdMag = (displFieldBgd(:,1).^2+displFieldBgd(:,2).^2).^0.5;
    if isempty(displFieldMag)
        dispDetec = 0;
    else
        % sort and take top 10 mags
        displFieldBgdMagsorted = sort(displFieldBgdMag);
        dispDetec = mean(displFieldMag)/mean(displFieldBgdMagsorted(1:floor(0.1*length(displFieldBgdMagsorted))));%1:10));%f
    end
end

% Load the forcefield
forceField=MD.getPackage(iPack).getProcess(4).loadChannelOutput;

% finding force at mesh location
org_fx = zeros(size(forceField(1).pos(:,1)));
org_fy = zeros(size(forceField(1).pos(:,1)));
nmfPoints = length(forceField(1).pos(:,1));
for k=1:nmfPoints
    [~,indcol_closest_x] = min(abs(x_mat_u(1,:)-forceField(1).pos(k,1)),[],2);
    [~,indrow_closest_y] = min(abs(y_mat_u(:,1)-forceField(1).pos(k,2)),[],1);
    row_bottom = max(1,indrow_closest_y-3);
    row_top = min(size(x_mat_u,1),indrow_closest_y+3);
    col_bottom = max(1,indcol_closest_x-3);
    col_top = min(size(y_mat_u,2),indcol_closest_x+3);
    loc_xmat = x_mat_u(row_bottom:row_top,col_bottom:col_top);
    loc_ymat = y_mat_u(row_bottom:row_top,col_bottom:col_top);
    loc_fx = force_x(row_bottom:row_top,col_bottom:col_top);
    loc_fy = force_y(row_bottom:row_top,col_bottom:col_top);
    org_fx(k) = interp2(loc_xmat,loc_ymat,loc_fx,forceField(1).pos(k,1),forceField(1).pos(k,2));
    if isnan(org_fx(k))
        org_fx(k) = force_x(indrow_closest_y,indcol_closest_x);
    end
    org_fy(k) = interp2(loc_xmat,loc_ymat,loc_fy,forceField(1).pos(k,1),forceField(1).pos(k,2));
    if isnan(org_fy(k))
        org_fy(k) = force_y(indrow_closest_y,indcol_closest_x);
    end
end

% errors in force field
meanForceError=nansum(((org_fx-forceField(1).vec(:,1)).^2+(org_fy-forceField(1).vec(:,2)).^2).^.5);
% heatmap creation and saving - i'll do it later

% force peak ratio
maskForce = ((x_mat_u-100).^2+(y_mat_u-150).^2)<=d/2;
% forceForceIdx = maskVectors(forceField(1).pos(:,1),forceField(1).pos(:,2),maskForce);
% make  a an interpolated TF image and get the peak force because force
% mesh is sparse
[fMap,XI,YI]=generateHeatmapFromField(forceField);
%new mask with XI and YI
maskForceXIYI = ((XI-100).^2+(YI-150).^2)<=d/2;

% if isempty(forceForceIdx)
%     peakForceRatio = 0;
% else
    x_vec = reshape(x_mat_u,[],1);
    y_vec = reshape(y_mat_u,[],1);
    force_x_vec = reshape(force_x,[],1);
    force_y_vec = reshape(force_y,[],1);
%     forceFieldForce = forceField(1).vec(forceForceIdx,:);
%     forceFieldMag = (forceFieldForce(:,1).^2+forceFieldForce(:,2).^2).^0.5;
    fMapFiltered = fMap.*maskForceXIYI;
    forceFieldMag = fMapFiltered(fMapFiltered>0);
    orgFieldForceIdx = maskVectors(x_vec,y_vec,maskForce);
    orgFieldForceMag = (force_x_vec(orgFieldForceIdx).^2+force_y_vec(orgFieldForceIdx).^2).^0.5;
    
    backgroundIdx = maskVectors(forceField(1).pos(:,1),forceField(1).pos(:,2),~bwmorph(maskForce,'dilate',10));
    forceFieldBgd = forceField(1).vec(backgroundIdx,:);
    forceFieldBgdMag = (forceFieldBgd(:,1).^2+forceFieldBgd(:,2).^2).^0.5;
    if isempty(forceFieldMag)
        peakForceRatio = 0;
        forceDetec = 0;
    else
        peakForceRatio = mean(forceFieldMag)/mean(orgFieldForceMag);
        forceDetec = mean(forceFieldMag)/max(forceFieldBgdMag);
    end
return

%% input parameters to be replaced with function inputs
% f=2000; %Pa
% d=10;
% minCorLength = 21;
% dataPath='/files/.retain-snapshots.d7d-w0d/LCCB/fsm/harvard/analysis/Sangyoon/Bead-tracking/singleForceTesting';
% 
% testSingleForce(f,d,minCorLength,dataPath)
