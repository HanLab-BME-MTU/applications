% runTestForceAdhProximity runs testForceAdhProximity
%% Simulations - simple enough simulation - only d
nExp = 5;
f = 400; % Pa
r = 6; % pixel
dmax = 32; % pixel
detectedL2 = false(nExp,dmax/2);
detectedL1 = false(nExp,dmax/2);
fm1L2 = zeros(nExp,dmax/2);
fm2L2 = zeros(nExp,dmax/2);
fm1L1 = zeros(nExp,dmax/2);
fm2L1 = zeros(nExp,dmax/2);
fMapL2 = cell(nExp,dmax/2);
fMapL1 = cell(nExp,dmax/2);
bead_x = cell(nExp);
bead_y = cell(nExp);
for epm=1:nExp
    ii=0;
    for d=2:2:dmax % distance between adhesions
        ii=ii+1;
        if ii==1
            method = 'L2';
            dataPath=['/project/cellbiology/gdanuser/adhesion/Sangyoon/TFM simulations/Proximity/exp' num2str(epm) 'f' num2str(f) 'd' num2str(d) method];
            [detectedL2(epm,ii),fm1L2(epm,ii),fm2L2(epm,ii),fMapL2{epm,ii},bead_x{epm}, bead_y{epm}, Av{epm}] = ...
                testForceAdhProximity(d,f,f,r,r,method,dataPath);
            method = 'L1';
            dataPath=['/project/cellbiology/gdanuser/adhesion/Sangyoon/TFM simulations/Proximity/exp' num2str(epm) 'f' num2str(f) 'd' num2str(d) method];
            [detectedL1(epm,ii),fm1L1(epm,ii),fm2L1(epm,ii),fMapL1{epm,ii}] = ...
                testForceAdhProximity(d,f,f,r,r,method,dataPath,bead_x{epm}, bead_y{epm}, Av{epm});
        else
            method = 'L2';
            dataPath=['/project/cellbiology/gdanuser/adhesion/Sangyoon/TFM simulations/Proximity/exp' num2str(epm) 'f' num2str(f) 'd' num2str(d) method];
            [detectedL2(epm,ii),fm1L2(epm,ii),fm2L2(epm,ii),fMapL2{epm,ii}] = ...
                testForceAdhProximity(d,f,f,r,r,method,dataPath,bead_x{epm}, bead_y{epm}, Av{epm});
            method = 'L1';
            dataPath=['/project/cellbiology/gdanuser/adhesion/Sangyoon/TFM simulations/Proximity/exp' num2str(epm) 'f' num2str(f) 'd' num2str(d) method];
            [detectedL1(epm,ii),fm1L1(epm,ii),fm2L1(epm,ii),fMapL1{epm,ii}] = ...
                testForceAdhProximity(d,f,f,r,r,method,dataPath,bead_x{epm}, bead_y{epm}, Av{epm});
        end
    end
end
%% save and plotting
save('/project/cellbiology/gdanuser/adhesion/Sangyoon/TFM simulations/Proximity/d_effect.mat')
% dataPath = '/hms/scratch1/sh268/singleForceTesting/f_vs_d/forceDetec';