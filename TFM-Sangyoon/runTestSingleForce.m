% runTestSingleForce tests effect of force magnitude and area of force
% application on correctness of displacement tracking and force
% reconstruction.
%% Simulations
d_err = zeros(20,10,3);
f_err = zeros(20,10,3);
ii=0;
jj=0;
kk=0;
p=0;
for f=200:200:4000 %Pa
    ii=ii+1;
    jj=0;
    for d=1:10
        jj=jj+1;
        kk=0;
        for cL = [11 21 41]
            p=p+1;
            kk=kk+1;
            dataPath=['/files/.retain-snapshots.d7d-w0d/LCCB/fsm/harvard/analysis/Sangyoon/Bead-tracking/singleForceTesting/f' num2str(f) 'd' num2str(d) 'cL' num2str(cL)];
            if p==1
                [d_err(ii,jj,kk),f_err(ii,jj,kk),bead_x, bead_y, Av]=testSingleForce(f,d,cL,dataPath); 
            else
                [d_err(ii,jj,kk),f_err(ii,jj,kk),~,~,~]=testSingleForce(f,d,cL,dataPath,bead_x, bead_y, Av);
            end
        end
    end
end
%% save
save('/files/.retain-snapshots.d7d-w0d/LCCB/fsm/harvard/analysis/Sangyoon/Bead-tracking/singleForceTesting/data.mat')
%%  visualize
[x,y] = meshgrid(10:-1:1,200:200:4000);
% figure,surf(x,y,newd_err(:,:,2))
% make interpolated surface
d_err_interp   = csapi({200:200:4000,10:-1:1},newd_err(:,:,2));
figure, fnplt( d_err_interp )
fu = 200:4000;
du = 1:.1:10;
[Cu,hu]=contour(fu,du,fnval(d_err_interp,{fu,du}).',20);
view(-90,90)
set(gca,'ydir','reverse');
set(hu,'ShowText','on')
% set(hu,'TextStep',get(hu,'LevelStep')*2)
set(gca,'Font','Arial');

colormap jet;
caxis(gca,'cdir','reverse')

%% force error
% figure,surf(x,y,newf_err(:,:,2))
f_err_interp   = csapi({200:200:4000,10:-1:1},newf_err(:,:,2));
% figure, fnplt( f_err_interp )
fu = 200:4000;
du = 1:.1:10;
figure,[Cf,hf]=contour(fu,du,fnval(f_err_interp,{fu,du}).',10);
view(-90,90)
set(gca,'ydir','reverse');
colormap jet;

%% original displacementfield when d is small
f=600; d=2; cL=21;
dataPath=['/files/.retain-snapshots.d7d-w0d/LCCB/fsm/harvard/analysis/Sangyoon/Bead-tracking/singleForceTesting/f' num2str(f) 'd' num2str(d) 'cL' num2str(cL)];
testSingleForce(f,d,cL,dataPath); 
% get the original displacement field
resultPath = [dataPath filesep 'Original' filesep 'data.mat'];
load(resultPath,'x_mat_u','y_mat_u','ux','uy')
generateHeatmapFromGridData(x_mat_u,y_mat_u,ux,uy,dataPath)
%% measured displacementfield when d is small
% get the measured displacement field
displPath = [dataPath filesep 'TFMPackage/displacementField'];
displFile = [dataPath filesep 'TFMPackage/displacementField/displField.mat'];
load(displFile)
generateHeatmapFromField(displField,displPath,.155)


%% original displacementfield when d is large
f=600; d=8; cL=21;
dataPath=['/files/.retain-snapshots.d7d-w0d/LCCB/fsm/harvard/analysis/Sangyoon/Bead-tracking/singleForceTesting/f' num2str(f) 'd' num2str(d) 'cL' num2str(cL)];
testSingleForce(f,d,cL,dataPath); 
%% get the original displacement field
resultPath = [dataPath filesep 'Original' filesep 'data.mat'];
load(resultPath,'x_mat_u','y_mat_u','ux','uy')
generateHeatmapFromGridData(x_mat_u,y_mat_u,ux,uy,dataPath)
%% measured displacementfield when d is small
% get the measured displacement field
displPath = [dataPath filesep 'TFMPackage/displacementField'];
displFile = [dataPath filesep 'TFMPackage/displacementField/displField.mat'];
load(displFile)
generateHeatmapFromField(displField,displPath,.71)

%% Bead tracking test with identical image stack
f=0; d=10; cL=21;
dataPath=['/files/.retain-snapshots.d7d-w0d/LCCB/fsm/harvard/analysis/Sangyoon/Bead-tracking/singleForceTesting/f' num2str(f) 'd' num2str(d) 'cL' num2str(cL)];
dataPath=['/Users/joshua2/Documents/PostdocResearch/Traction Force/corrTrackContWind/f' num2str(f) 'd' num2str(d) 'cL' num2str(cL)];
[derr0f, ferr0f] = testSingleForce(f,d,cL,dataPath); 
%% measured displacementfield for zero force
% get the measured displacement field
displPath = [dataPath filesep 'TFMPackage/displacementField'];
displFile = [dataPath filesep 'TFMPackage/displacementField/displField.mat'];
load(displFile)
generateHeatmapFromField(displField,displPath,.15)
