% runTestSingleForce tests effect of force magnitude and area of force
% application on correctness of displacement tracking and force
% reconstruction.
%% Simulations - initialization  for f and d
nExp = 5;
d_err_FDAdhL1new = zeros(20,10,nExp);
d_err_FDBGL1new = zeros(20,10,nExp);
f_err_FDADhL1new = zeros(20,10,nExp);
f_err_FDBGL1new = zeros(20,10,nExp);
dispDetec_FDL1new = zeros(20,10,nExp);
forceDetec_FDL1new = zeros(20,10,nExp);
pFR_FDL1new = zeros(20,10,nExp);
beadsOnAdhnew = zeros(nExp,1);
cL = 15;%[9 15 21]
% kk=0;
% simulation for f and d (L1 with 10% noise, old tracking)
for epm=1:nExp
    p=0;
    ii=0;
    for f=200:200:4000 %Pa
        ii=ii+1;
        jj=0;
        for d=2:2:20
            jj=jj+1;
    %         kk=0;
            p=p+1;
%             kk=kk+1;
            dataPath=['/hms/scratch1/sh268/singleForceTesting/f_vs_d/simulations/exp' num2str(epm) 'f' num2str(f) 'd' num2str(d) 'cL' num2str(cL) 'L1'];
            if p==1
                [d_err_FDAdhL1new(ii,jj,epm),d_err_FDBGL1new(ii,jj,epm),dispDetec_FDL1new(ii,jj,epm),...
                    f_err_FDADhL1new(ii,jj,epm),f_err_FDBGL1new(ii,jj,epm),pFR_FDL1new(ii,jj,epm),forceDetec_FDL1new(ii,jj,epm),...
                    beadsOnAdhnew(epm),bead_xL2{epm}, bead_yL2{epm}, AvL2{epm}]= ...
                    testSingleForce(f,d,cL,dataPath,[],[],[],'1NormReg');
            else
                [d_err_FDAdhL1new(ii,jj,epm),d_err_FDBGL1new(ii,jj,epm),dispDetec_FDL1new(ii,jj,epm),...
                    f_err_FDADhL1new(ii,jj,epm),f_err_FDBGL1new(ii,jj,epm),pFR_FDL1new(ii,jj,epm),forceDetec_FDL1new(ii,jj,epm),...
                    beadsOnAdhnew(epm)]= ...
                    testSingleForce(f,d,cL,dataPath,bead_xL2{epm}, bead_yL2{epm}, AvL2{epm},'1NormReg');
            end
%                 dataPath=['/files/.retain-snapshots.d7d-w0d/LCCB/fsm/harvard/analysis/Sangyoon/Bead-tracking/singleForceTesting/f_vs_d/simulations/exp' num2str(epm) 'f' num2str(f) 'd' num2str(d) 'cL' num2str(cL)];
%                 if p==1
%                     [d_err(ii,jj,epm),dispDetec(ii,jj,epm),f_err(ii,jj,epm),peakForceRatio(ii,jj,epm),forceDetec(ii,jj,epm),bead_x, bead_y, Av]=testSingleForce(f,d,cL,dataPath); 
%     %                 [d_err(ii,jj,kk),dispDetec(ii,jj,kk),f_err(ii,jj,kk),peakForceRatio(ii,jj,kk),bead_x, bead_y, Av]=testSingleForce(f,d,cL,dataPath); 
%                 else
%                     [d_err(ii,jj,epm),dispDetec(ii,jj,epm),f_err(ii,jj,epm),peakForceRatio(ii,jj,epm),forceDetec(ii,jj,epm),~,~,~]=testSingleForce(f,d,cL,dataPath,bead_x, bead_y, Av);
%     %                 [d_err(ii,jj,kk),dispDetec(ii,jj,kk),f_err(ii,jj,kk),peakForceRatio(ii,jj,kk),~,~,~]=testSingleForce(f,d,cL,dataPath,bead_x, bead_y, Av);
%                 end
        end
    end
end
%% save and plotting
save('/hms/scratch1/sh268/singleForceTesting/f_vs_d/FvsD_L1.mat')
% dataPath = '/hms/scratch1/sh268/singleForceTesting/f_vs_d/forceDetec';
%% showing for L1 heatmap
f=200:200:4000;
d=2:2:12;
dataPath = '/Users/joshua2/Documents/PostdocResearch/Traction Force/Manuscript/Draft1.2/Data/SingleForceTesting/forceDetec_L1';
visualizeError(f,d,forceDetec_FDL1new(:,1:6,:),dataPath,'contourf_with_level1_2',3)
%% Simulations - initialization  for f and d
d_err_FDAdhL2new = zeros(20,10,nExp);
d_err_FDBGL2new = zeros(20,10,nExp);
f_err_FDADhL2new = zeros(20,10,nExp);
f_err_FDBGL2new = zeros(20,10,nExp);
dispDetec_FDL2new = zeros(20,10,nExp);
forceDetec_FDL2new = zeros(20,10,nExp);
pFR_FDL2new = zeros(20,10,nExp);
cL = 15;%[9 15 21]
% kk=0;
% simulation for f and d (L2 with 10% noise, old tracking)
for epm=1:nExp
    p=0;
    ii=0;
    for f=200:200:4000 %Pa
        ii=ii+1;
        jj=0;
        for d=2:2:20
            jj=jj+1;
            p=p+1;
            dataPath=['/hms/scratch1/sh268/singleForceTesting/f_vs_d/simulations/exp' num2str(epm) 'f' num2str(f) 'd' num2str(d) 'cL' num2str(cL) 'L1'];
            storagePath=['/hms/scratch1/sh268/singleForceTesting/f_vs_d/simulations/exp' num2str(epm) 'f' num2str(f) 'd' num2str(d) 'cL' num2str(cL) 'L2new'];
            [d_err_FDAdhL2new(ii,jj,epm),d_err_FDBGL2new(ii,jj,epm),dispDetec_FDL2new(ii,jj,epm),...
                    f_err_FDADhL2new(ii,jj,epm),f_err_FDBGL2new(ii,jj,epm),pFR_FDL2new(ii,jj,epm),forceDetec_FDL2new(ii,jj,epm),...
                    beadsOnAdhnew(epm)]= testSingleForceChange(d,dataPath,storagePath,'QR',1e-6);
        end
    end
end
%% save
clear d_err_FDAdhL1new d_err_FDBGL1new f_err_FDADhL1new f_err_FDBGL1new dispDetec_FDL1new forceDetec_FDL1new pFR_FDL
save('/hms/scratch1/sh268/singleForceTesting/f_vs_d/FvsD_L2new.mat')
%% reanalyze for L2
d_err_FDAdhL2new = zeros(20,10,nExp);
d_err_FDBGL2new = zeros(20,10,nExp);
f_err_FDADhL2new = zeros(20,10,nExp);
f_err_FDBGL2new = zeros(20,10,nExp);
dispDetec_FDL2new = zeros(20,10,nExp);
forceDetec_FDL2new = zeros(20,10,nExp);
pFR_FDL2new = zeros(20,10,nExp);
cL = 15;%[9 15 21]
L=length(num2str(nExp));
strg=sprintf('%%.%dd',L);
backSpc =repmat('\b',1,L);
% kk=0;
% simulation for f and d (L2 with 10% noise, old tracking)
for epm=1:nExp
    fprintf(1,[strg ' ...'],epm);
    p=0;
    ii=0;
    for f=200:200:4000 %Pa
        ii=ii+1;
        jj=0;
        for d=2:2:20
            jj=jj+1;
            p=p+1;
            dataPath=['/hms/scratch1/sh268/singleForceTesting/f_vs_d/simulations/exp' num2str(epm) 'f' num2str(f) 'd' num2str(d) 'cL' num2str(cL) 'L1'];
            storagePath=['/hms/scratch1/sh268/singleForceTesting/f_vs_d/simulations/exp' num2str(epm) 'f' num2str(f) 'd' num2str(d) 'cL' num2str(cL) 'L2new'];
            [d_err_FDAdhL2new(ii,jj,epm),d_err_FDBGL2new(ii,jj,epm),dispDetec_FDL2new(ii,jj,epm),...
                    f_err_FDADhL2new(ii,jj,epm),f_err_FDBGL2new(ii,jj,epm),pFR_FDL2new(ii,jj,epm),forceDetec_FDL2new(ii,jj,epm),...
                    beadsOnAdhnew(epm)]= analyzeSingleForceDataChange(d,dataPath,storagePath);
        end
    end
    fprintf(1,[backSpc '\b\b\b\b']);
end
%% save
clear d_err_FDAdhL1new d_err_FDBGL1new f_err_FDADhL1new f_err_FDBGL1new dispDetec_FDL1new forceDetec_FDL1new pFR_FDL
save('/hms/scratch1/sh268/singleForceTesting/f_vs_d/FvsD_L2new.mat')
%% Simulations - with pixel-wise tracking
d_err_FDAdhL2old = zeros(20,10,nExp);
d_err_FDBGL2old = zeros(20,10,nExp);
f_err_FDADhL2old = zeros(20,10,nExp);
f_err_FDBGL2old = zeros(20,10,nExp);
dispDetec_FDL2old = zeros(20,10,nExp);
forceDetec_FDL2old = zeros(20,10,nExp);
pFR_FDL2old = zeros(20,10,nExp);
cL = 15;%[9 15 21]
% kk=0;
% simulation for f and d (L2 with 10% noise, old tracking)
for epm=1:nExp
    p=0;
    ii=0;
    for f=200:200:4000 %Pa
        ii=ii+1;
        jj=0;
        for d=2:2:20
            jj=jj+1;
            p=p+1;
            dataPath=['/hms/scratch1/sh268/singleForceTesting/f_vs_d/simulations/exp' num2str(epm) 'f' num2str(f) 'd' num2str(d) 'cL' num2str(cL) 'L2old'];
            [d_err_FDAdhL2old(ii,jj,epm),d_err_FDBGL2old(ii,jj,epm),dispDetec_FDL2old(ii,jj,epm),...
                f_err_FDADhL2old(ii,jj,epm),f_err_FDBGL2old(ii,jj,epm),pFR_FDL2old(ii,jj,epm),forceDetec_FDL2old(ii,jj,epm),...
                beadsOnAdhnew(epm)]= ...
                testSingleForce(f,d,cL,dataPath,bead_xL2{epm}, bead_yL2{epm}, AvL2{epm},'QR');
        end
    end
end
save('/hms/scratch1/sh268/singleForceTesting/f_vs_d/FvsD_L2old.mat')
    %% show heatmap
f=200:200:4000;
d=2:2:12;
% dataPath = '/hms/scratch1/sh268/singleForceTesting/f_vs_d/forceDetec_L2new';
dataPath = '/Users/joshua2/Documents/PostdocResearch/Traction Force/Manuscript/Draft1.2/Data/SingleForceTesting/forceDetec_L2new';
visualizeError(f,d,forceDetec_FDL2new(:,1:6,:),dataPath,'contourf_with_level1_2',3)
%% for L2 old
f=200:200:4000;
d=2:2:12;
% dataPath = '/hms/scratch1/sh268/singleForceTesting/f_vs_d/forceDetec_L2new';
dataPath = '/Users/joshua2/Documents/PostdocResearch/Traction Force/Manuscript/Draft1.2/Data/SingleForceTesting/forceDetec_L2old';
visualizeError(f,d,forceDetec_FDL2old(:,1:6,:),dataPath,'contourf_with_level1_2',3)

%% load
load('/files/.retain-snapshots.d7d-w0d/LCCB/fsm/harvard/analysis/Sangyoon/Bead-tracking/singleForceTesting/data.mat')
%% reanalyze for dispDetec
for epm=1%:nExp
    p=0;
    ii=0;
    for f=200:200:4000 %Pa
        ii=ii+1;
        jj=0;
        for d=2:2:20
            jj=jj+1;
    %         kk=0;
            for cL = 15%[9 15 21]
                p=p+1;
    %             kk=kk+1;
                dataPath=['/files/.retain-snapshots.d7d-w0d/LCCB/fsm/harvard/analysis/Sangyoon/Bead-tracking/singleForceTesting/f_vs_d/simulations/exp' num2str(epm) 'f' num2str(f) 'd' num2str(d) 'cL' num2str(cL)];
                [~,dispDetec(ii,jj,epm),~,~,~]=analyzeSingleForceData(f,d,cL,dataPath);
            end
        end
    end
    disp(['epm=' num2str(epm)])
end
%%  visualize (contour)
% load('/files/.retain-snapshots.d7d-w0d/LCCB/fsm/harvard/analysis/Sangyoon/Bead-tracking/singleForceTesting_old/data.mat')
f=200:200:4000;
d=2:2:20;
% dispDetec(:,:,1) = meshgrid(d,f);
% dataPath = '/files/.retain-snapshots.d7d-w0d/LCCB/fsm/harvard/analysis/Sangyoon/Bead-tracking/singleForceTesting/f_vs_d/dispDetec';
% visualizeError(f,d,dispDetec,dataPath)
dataPath = '/hms/scratch1/sh268/singleForceTesting/f_vs_d/forceDetec';
visualizeError(f,d,forceDetec_FDL1new,dataPath,'pcolor_with_level1line',100)

dataPath = '/hms/scratch1/sh268/singleForceTesting/f_vs_d/peakForceRatio';
visualizeError(f,d,pFR_FDL1new,dataPath,'contourf')
%% for old trackstackflow
nExp = 10;
d_err_old = zeros(10,nExp);
f_err_old = zeros(10,nExp);
dispDetec_old = zeros(10,nExp);
forceDetec_old = zeros(10,nExp);
peakForceRatio_old = zeros(10,nExp);
beadsOnAdhold = zeros(nExp,1);
f=1000;
cL = 15;
for epm=1:nExp
    p=0;
    jj=0;
    tstart = cputime;
    for d=2:2:20
        jj=jj+1;
        p=p+1;
%         dataPath=['/files/.retain-snapshots.d7d-w0d/LCCB/fsm/harvard/analysis/Sangyoon/Bead-tracking/singleForceTesting/f_vs_d/simulations/exp' num2str(epm) 'f' num2str(f) 'd' num2str(d) 'cL' num2str(cL) 'interp'];
        dataPath=['/home/sh268/files/LCCB/fsm/harvard/analysis/Sangyoon/Bead-tracking/singleForceTesting/f_vs_d/simulations/exp' num2str(epm) 'f' num2str(f) 'd' num2str(d) 'cL' num2str(cL) 'interp'];
        if p==1
            [d_err_old(jj,epm),dispDetec_old(jj,epm),f_err_old(jj,epm),peakForceRatio_old(jj,epm),forceDetec_old(jj,epm),beadsOnAdhold(epm),bead_x{epm}, bead_y{epm}, Av{epm}]= testSingleForce(f,d,cL,dataPath,[],[],[],'QR');
        else
            [d_err_old(jj,epm),dispDetec_old(jj,epm),f_err_old(jj,epm),peakForceRatio_old(jj,epm),forceDetec_old(jj,epm),~,~,~,~]=testSingleForce(f,d,cL,dataPath,bead_x{epm}, bead_y{epm}, Av{epm},'QR');
        end
    end
    disp(['Experiment ' num2str(epm) ' is done! ' num2str(cputime-tstart) ' is elapsed. ' num2str((cputime-tstart)*(nExp-epm)) ' sec is expected.'])
end
%% for new trackstackflow
% change setting for trackStackFlow
d_err_new = zeros(10,nExp);
f_err_new = zeros(10,nExp);
dispDetec_new = zeros(10,nExp);
forceDetec_new = zeros(10,nExp);
peakForceRatio_new = zeros(10,nExp);
beadsOnAdhnew = zeros(nExp,1);
f=1000;
cL = 15;
for epm=1:nExp
    p=0;
    jj=0;
    for d=2:2:20
        jj=jj+1;
        p=p+1;
%         dataPath=['/files/.retain-snapshots.d7d-w0d/LCCB/fsm/harvard/analysis/Sangyoon/Bead-tracking/singleForceTesting/f_vs_d/simulations/exp' num2str(epm) 'f' num2str(f) 'd' num2str(d) 'cL' num2str(cL) 'interp'];
        dataPath=['/home/sh268/files/LCCB/fsm/harvard/analysis/Sangyoon/Bead-tracking/singleForceTesting/f_vs_d/simulations/exp' num2str(epm) 'f' num2str(f) 'd' num2str(d) 'cL' num2str(cL)];
        [d_err_new(jj,epm),dispDetec_new(jj,epm),f_err_new(jj,epm),peakForceRatio_new(jj,epm),forceDetec_new(jj,epm),beadsOnAdhnew(epm),~,~,~]=testSingleForce(f,d,cL,dataPath,bead_x{epm}, bead_y{epm}, Av{epm},'QR');
    end
end
%% save the data
save('/home/sh268/files/LCCB/fsm/harvard/analysis/Sangyoon/Bead-tracking/singleForceTesting/f_vs_d/forceDetecVsDF1000.mat')
%% d vs forceDetec for f=1000Pa
d=2:2:20;
figure, plot(d,mean(forceDetec(5,:,:),3)), hold on, plot(d,mean(forceDetec_old,2),'r')
figure, plot(d,mean(f_err(5,:,:),3)), hold on, plot(d,mean(f_err_old,2),'r')
figure, plot(d,mean(peakForceRatio(5,:,:),3)), hold on, plot(d,mean(peakForceRatio_old,2),'r')
%% d vs forceDetec for f=1000Pa
d=2:2:20;
figure, plot(d,mean(forceDetec_FDL2old(5,:,:),3))
figure, plot(d,mean(forceDetec_FDL2old(2,:,:),3))
%% surf
figure, surf(x,y,mean(peakForceRatio(:,:,1:2),3)), title('peak force ratio')
figure, surf(x,y,mean(f_err(:,:,1:2),3)), title('force RMS error')
figure, surf(x,y,mean(d_err(:,:,1:2),3)), title('displacement RMS error')
figure, surf(x,y,mean(dispDetec(:,:,1:2),3)), title('detectability')

%% template size on force error, using L2 0th with a new method
nExp = 5;
f=2000;
d=6;
d_err_tsQRnew = zeros(19,nExp);
f_err_tsQRnew = zeros(19,nExp);
dispDetec_tsQRnew = zeros(19,nExp);
forceDetec_tsQRnew = zeros(19,nExp);
peakForceRatio_tsQRnew = zeros(19,nExp);
for epm=1:nExp
    p=0;
    for cL = 7:2:43
        p=p+1;
        dataPath=['/home/sh268/files/LCCB/fsm/harvard/analysis/Sangyoon/Bead-tracking/singleForceTesting/f_vs_d/simulations/exp' num2str(epm) 'f' num2str(f) 'd' num2str(d) 'cL' num2str(cL)];
        [d_err_tsQRnew(p,epm),dispDetec_tsQRnew(p,epm),f_err_tsQRnew(p,epm),peakForceRatio_tsQRnew(p,epm),forceDetec_tsQRnew(p,epm),~,~,~,~]=testSingleForce(f,d,cL,dataPath,bead_x{epm}, bead_y{epm}, Av{epm},'QR');
    end
end
%% template size on force error, using L2 0th with a old method
nExp = 5;
f=2000;
d=6;
d_err_tsQRold = zeros(19,nExp);
f_err_tsQRold = zeros(19,nExp);
dispDetec_tsQRold = zeros(19,nExp);
forceDetec_tsQRold = zeros(19,nExp);
peakForceRatio_tsQRold = zeros(19,nExp);
for epm=1:nExp
    p=0;
    for cL = 7:2:43
        p=p+1;
        dataPath=['/home/sh268/files/LCCB/fsm/harvard/analysis/Sangyoon/Bead-tracking/singleForceTesting/f_vs_d/simulations/exp' num2str(epm) 'f' num2str(f) 'd' num2str(d) 'cL' num2str(cL) ' interp'];
        [d_err_tsQRold(p,epm),dispDetec_tsQRold(p,epm),f_err_tsQRold(p,epm),peakForceRatio_tsQRold(p,epm),forceDetec_tsQRold(p,epm),~,~,~,~]=testSingleForce(f,d,cL,dataPath,bead_x{epm}, bead_y{epm}, Av{epm},'QR');
    end
end
%% reanalyze for cL effect with new tracking
nExp = 5;
f=2000;
d=6;
d_err_cLAdhnew1 = zeros(19,nExp);
d_err_cLBGnew1 = zeros(19,nExp);
f_err_cLAdhnew1 = zeros(19,nExp);
f_err_cLBGnew1 = zeros(19,nExp);
dispDetec_tsQRnew1 = zeros(19,nExp);
forceDetec_tsQRnew1 = zeros(19,nExp);
peakForceRatio_tsQRnew1 = zeros(19,nExp);
for epm=1:nExp
    p=0;
    for cL = 7:2:43
        p=p+1;
        dataPath=['/home/sh268/files/LCCB/fsm/harvard/analysis/Sangyoon/Bead-tracking/singleForceTesting/f_vs_d/simulations/exp' num2str(epm) 'f' num2str(f) 'd' num2str(d) 'cL' num2str(cL)];
        [d_err_cLAdhnew1(p,epm),d_err_cLBGnew1(p,epm),dispDetec_tsQRnew1(p,epm),f_err_cLAdhnew1(p,epm),f_err_cLBGnew1(p,epm),peakForceRatio_tsQRnew1(p,epm),forceDetec_tsQRnew1(p,epm),~]=analyzeSingleForceData(d,dataPath);
    end
end
%% just for peakForceRatio
[d_err_atemp,d_err_btemp,dispDetec_temp,f_err_atemp,f_err_btemp,peakForceRatio_temp,forceDetec_temp,~]=analyzeSingleForceData(d,dataPath);
%% reanalyze for cL effect with old tracking
nExp = 5;
f=2000;
d=6;
d_err_cLAdhold1 = zeros(19,nExp);
d_err_cLBGold1 = zeros(19,nExp);
f_err_cLAdhold1 = zeros(19,nExp);
f_err_cLBGold1 = zeros(19,nExp);
dispDetec_tsQRold1 = zeros(19,nExp);
forceDetec_tsQRold1 = zeros(19,nExp);
peakForceRatio_tsQRold1 = zeros(19,nExp);
for epm=1:nExp
    p=0;
    for cL = 7:2:43
        p=p+1;
        dataPath=['/home/sh268/files/LCCB/fsm/harvard/analysis/Sangyoon/Bead-tracking/singleForceTesting/f_vs_d/simulations/exp' num2str(epm) 'f' num2str(f) 'd' num2str(d) 'cL' num2str(cL) ' interp'];
        [d_err_cLAdhold1(p,epm),d_err_cLBGold1(p,epm),dispDetec_tsQRold1(p,epm),f_err_cLAdhold1(p,epm),f_err_cLBGold1(p,epm),peakForceRatio_tsQRold1(p,epm),forceDetec_tsQRold1(p,epm),~]=analyzeSingleForceData(d,dataPath);
    end
end
%% save the data
save('/home/sh268/files/LCCB/fsm/harvard/analysis/Sangyoon/Bead-tracking/singleForceTesting/f_vs_d/cLeffect.mat')
%% template size on force error with 10% noise + less gap, using L2 0th with a old method
nExp = 5;
f=2000;
d=6;
d_err_cLAdhold10 = zeros(19,nExp);
d_err_cLBGold10 = zeros(19,nExp);
f_err_cLAdhold10 = zeros(19,nExp);
f_err_cLBGold10 = zeros(19,nExp);
dispDetec_tsQRold10 = zeros(19,nExp);
forceDetec_tsQRold10 = zeros(19,nExp);
peakForceRatio_tsQRold10 = zeros(19,nExp);
for epm=1:nExp
    p=0;
    for cL = 7:2:43
        p=p+1;
        dataPath=['/home/sh268/files/LCCB/fsm/harvard/analysis/Sangyoon/Bead-tracking/singleForceTesting/f_vs_d/simulations/exp' num2str(epm) 'f' num2str(f) 'd' num2str(d) 'cL' num2str(cL) 'n10 interp'];
        if p==1
            [d_err_cLAdhold10(p,epm),d_err_cLBGold10(p,epm),dispDetec_tsQRold10(p,epm),f_err_cLAdhold10(p,epm),f_err_cLBGold10(p,epm),peakForceRatio_tsQRold10(p,epm),forceDetec_tsQRold10(p,epm),beadsOnAdhold(epm),bead_x10{epm}, bead_y10{epm}, Av10{epm}]= ...
                testSingleForce(f,d,cL,dataPath,[],[],[],'QR');
        else
            [d_err_cLAdhold10(p,epm),d_err_cLBGold10(p,epm),dispDetec_tsQRold10(p,epm),f_err_cLAdhold10(p,epm),f_err_cLBGold10(p,epm),peakForceRatio_tsQRold10(p,epm),forceDetec_tsQRold10(p,epm),~,~,~,~]...
                =testSingleForce(f,d,cL,dataPath,bead_x10{epm}, bead_y10{epm}, Av10{epm},'QR');
        end
    end
end
%% template size on force error with 10% noise + less gap, using L2 0th with a new method
nExp = 5;
f=2000;
d=6;
d_err_cLAdhnew10 = zeros(19,nExp);
d_err_cLBGnew10 = zeros(19,nExp);
f_err_cLAdhnew10 = zeros(19,nExp);
f_err_cLBGnew10 = zeros(19,nExp);
dispDetec_tsQRnew10 = zeros(19,nExp);
forceDetec_tsQRnew10 = zeros(19,nExp);
peakForceRatio_tsQRnew10 = zeros(19,nExp);
for epm=1:nExp
    p=0;
    for cL = 7:2:43
        p=p+1;
        dataPath=['/home/sh268/files/LCCB/fsm/harvard/analysis/Sangyoon/Bead-tracking/singleForceTesting/f_vs_d/simulations/exp' num2str(epm) 'f' num2str(f) 'd' num2str(d) 'cL' num2str(cL) 'n10'];

        [d_err_cLAdhnew10(p,epm),d_err_cLBGnew10(p,epm),dispDetec_tsQRnew10(p,epm),f_err_cLAdhnew10(p,epm),f_err_cLBGnew10(p,epm),peakForceRatio_tsQRnew10(p,epm),forceDetec_tsQRnew10(p,epm),~,~,~,~]=...
            testSingleForce(f,d,cL,dataPath,bead_x10{epm}, bead_y10{epm}, Av10{epm},'QR');
    end
end
%% plotting cL vs. force error
cL= 9:2:43;
figure, plot(cL,f_err_tsQRnew)
figure, plot(cL,peakForceRatio_tsQRnew)
figure, plot(cL,forceDetec_tsQRnew)
figure, plot(cL,d_err_tsQRnew)
hold on
plot(cL,mean(d_err_tsQRnew,2),'k','LineWidth',2)
%% showing displacement maps for each condition
for cL = 9:2:41
    dataPath=['/files/.retain-snapshots.d7d-w0d/LCCB/fsm/harvard/analysis/Sangyoon/Bead-tracking/singleForceTesting/cL_effect/simulations/f' num2str(f) 'd' num2str(d) 'cL' num2str(cL)];
    displPath = [dataPath filesep 'TFMPackage/displacementField'];
    displFile = [dataPath filesep 'TFMPackage/displacementField/displField.mat'];
    load(displFile)
    generateHeatmapFromField(displField,displPath,3.7);
end
%% get the original displacement field
resultPath = [dataPath filesep 'Original' filesep 'data.mat'];
load(resultPath,'x_mat_u','y_mat_u','ux','uy')
generateHeatmapFromGridData(x_mat_u,y_mat_u,ux,uy,dataPath)

%% showing force maps for each condition
for cL = 9:2:41
    dataPath=['/files/.retain-snapshots.d7d-w0d/LCCB/fsm/harvard/analysis/Sangyoon/Bead-tracking/singleForceTesting/cL_effect/simulations/f' num2str(f) 'd' num2str(d) 'cL' num2str(cL)];
    generateHeatmapFromTFMPackage(dataPath,16,false,2000)
end
resultPath = [dataPath filesep 'Original' filesep 'data.mat'];
load(resultPath,'x_mat_u','y_mat_u','force_x','force_y')
generateHeatmapFromGridData(x_mat_u,y_mat_u,force_x,force_y,dataPath)

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
epm=2; f=600; d=4; cL=15;
% dataPath=['/files/.retain-snapshots.d7d-w0d/LCCB/fsm/harvard/analysis/Sangyoon/Bead-tracking/singleForceTesting/f_vs_d/simulations/exp' num2str(epm) 'f' num2str(f) 'd' num2str(d) 'cL' num2str(cL) 'interp'];
% dataPath=['/files/.retain-snapshots.d7d-w0d/LCCB/fsm/harvard/analysis/Sangyoon/Bead-tracking/singleForceTesting/f_vs_d/simulations/exp' num2str(epm) 'f' num2str(f) 'd' num2str(d) 'cL' num2str(cL)];
dataPath=['/home/sh268/files/LCCB/fsm/harvard/analysis/Sangyoon/Bead-tracking/singleForceTesting/f_vs_d/simulations/exp' num2str(epm) 'f' num2str(f) 'd' num2str(d) 'cL' num2str(cL)];
[d_err_small,dispDetec_small,f_err_small,peakForceRatio_small,forceDetec_small,bead_x, bead_y, Av]= testSingleForce(f,d,cL,dataPath,[],[],[],'QR');
            
%% get the original displacement field
resultPath = [dataPath filesep 'Original' filesep 'data.mat'];
load(resultPath,'x_mat_u','y_mat_u','ux','uy')
generateHeatmapFromGridData(x_mat_u,y_mat_u,ux,uy,dataPath,100)
%% show original force field
load(resultPath,'force_x','force_y')
generateHeatmapFromGridData(x_mat_u,y_mat_u,force_x,force_y,[dataPath '/Original forcefield'],1000,140,220)
%% measured displacementfield when d is small
% get the measured displacement field
displPath = [dataPath filesep 'TFMPackage/displacementField'];
displFile = [dataPath filesep 'TFMPackage/displacementField/displField.mat'];
load(displFile)
generateHeatmapFromField(displField,displPath,0.6,'uDefinedRYG',140,220);
% generateHeatmapFromField(displField,displPath,0.25,'cool');
%% measured forcemap when d is small
forcePath = [dataPath filesep 'TFMPackage/forceField'];
forceFile = [dataPath filesep 'TFMPackage/forceField/forceField.mat'];
load(forceFile)
generateHeatmapFromField(forceField,forcePath,1000,'jet',140,220);

%% original displacementfield when d is large
epm=2; f=1000; d=20; cL=15;
% dataPath=['/files/.retain-snapshots.d7d-w0d/LCCB/fsm/harvard/analysis/Sangyoon/Bead-tracking/singleForceTesting/f_vs_d/simulations/exp' num2str(epm) 'f' num2str(f) 'd' num2str(d) 'cL' num2str(cL) 'interp'];
dataPath=['/files/.retain-snapshots.d7d-w0d/LCCB/fsm/harvard/analysis/Sangyoon/Bead-tracking/singleForceTesting/f_vs_d/simulations/exp' num2str(epm) 'f' num2str(f) 'd' num2str(d) 'cL' num2str(cL)];
dataPath=['/home/sh268/files/LCCB/fsm/harvard/analysis/Sangyoon/Bead-tracking/singleForceTesting/f_vs_d/simulations/exp' num2str(epm) 'f' num2str(f) 'd' num2str(d) 'cL' num2str(cL) 'interp'];
[d_err_large,dispDetec_large,f_err_large,peakForceRatio_large,forceDetec_large]= testSingleForce(f,d,cL,dataPath,[],[],[],'QR');
%% show original force field
resultPath = [dataPath filesep 'Original' filesep 'data.mat'];
load(resultPath,'x_mat_u','y_mat_u','ux','uy','force_x','force_y')
generateHeatmapFromGridData(x_mat_u,y_mat_u,force_x,force_y,[dataPath '/Original forcefield'])
%% get the original displacement field
resultPath = [dataPath filesep 'Original' filesep 'data.mat'];
load(resultPath,'x_mat_u','y_mat_u','ux','uy')
generateHeatmapFromGridData(x_mat_u,y_mat_u,ux,uy,dataPath)
%% measured displacementfield when d is large
% get the measured displacement field
displPath = [dataPath filesep 'TFMPackage/displacementField'];
displFile = [dataPath filesep 'TFMPackage/displacementField/displField.mat'];
load(displFile)
generateHeatmapFromField(displField,displPath,2.8,'uDefinedRYG',140,220);
% generateHeatmapFromField(displField,displPath,2.6);
%% measured forcemap when d is large
forcePath = [dataPath filesep 'TFMPackage/forceField'];
forceFile = [dataPath filesep 'TFMPackage/forceField/forceField.mat'];
load(forceFile)
generateHeatmapFromField(forceField,forcePath,1000,'jet',140,220);

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

%% L-curve analysis for /home/sh268/files/LCCB/fsm/harvard/analysis/Sangyoon/Bead-tracking/singleForceTesting/f_vs_d/simulations/exp3f1000d4cL7
dataPath = '/home/sh268/files/LCCB/fsm/harvard/analysis/Sangyoon/Bead-tracking/singleForceTesting/f_vs_d/simulations/exp1f1000d4cL15';

disp('calculating L-curve with L2 0th...')
load([dataPath '/TFMPackage/forceField/BEMParams.mat']);
MpM = M'*M;
[eyeWeights,~] =getGramMatrix(forceMesh);
% original force at force base nodes
xminf = forceMesh.basis(1).node(1);
yminf = forceMesh.basis(1).node(2);
gridSpacingf = forceMesh.basis(2).node(2)-forceMesh.basis(1).node(2);
xmaxf = forceMesh.basis(end).node(1);
ymaxf = forceMesh.basis(end).node(2);

force_x_f = force_x(yminf:gridSpacingf:ymaxf,xminf:gridSpacingf:xmaxf);
force_y_f = force_y(yminf:gridSpacingf:ymaxf,xminf:gridSpacingf:xmaxf);

force_x_vec_f=reshape(force_x_f,[],1);
force_y_vec_f=reshape(force_y_f,[],1);
force_0=vertcat(force_x_vec_f,force_y_vec_f);

alphas=10.^(-9:0.125:-3);
rho=zeros(length(alphas),1);
eta=zeros(length(alphas),1);
fErr=zeros(length(alphas),1);
msparse=zeros(size(M,2),length(alphas));
maxIter = 10;
tolx = 2e-2;
tolr = 1e-7;

for i=1:length(alphas);
  msparse(:,i)=iterativeL1Regularization(M,MpM,u,eyeWeights,alphas(i),maxIter,tolx,tolr);
  rho(i)=norm(M*msparse(:,i)-u);
  eta(i)=norm(msparse(:,i),1);
  % force error
  fErr(i)=norm(msparse(:,i)-force_0);
  disp([num2str(i) ' out of ' num2str(length(alphas))]);
end

%% Find the corner of the Tikhonov L-curve
[reg_corner,ireg_corner,kappa]=regParamSelecetionLcurve(rho,eta,alphas);
[~,fminIdx]=min(fErr);

save([dataPath '/LcurveL1-0th.mat'],'rho','eta','fErr','reg_corner','ireg_corner','alphas','msparse','fminIdx','-v7.3');
%% ireg_corner
[xgrid,ygrid]=meshgrid(xminf:gridSpacingf:xmaxf,yminf:gridSpacingf:ymaxf);
[fx,fy,x_out,y_out]=calcForcesFromCoef(forceMesh,msparse(:,ireg_corner),xgrid,ygrid,'new');
generateHeatmapFromGridData(x_out,y_out,fx,fy,[dataPath filesep 'L1forcemap at Lcorner'])
%% fminIdx
[fx,fy,x_out,y_out]=calcForcesFromCoef(forceMesh, msparse(:,fminIdx),xgrid,ygrid,'new');
generateHeatmapFromGridData(x_out,y_out,fx,fy,[dataPath filesep 'L1forcemap at fErr minimum'])
