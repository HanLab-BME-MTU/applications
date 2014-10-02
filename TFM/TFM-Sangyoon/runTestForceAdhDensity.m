% runTestForceAdhProximity runs testForceAdhProximity
%% Simulations - simple enough simulation - only d
xmax=200;
ymax=200;
nExp = 1;
f = 400; % Pa
r = 3; % radius in pixel

posx_min = 4*r; posx_max = xmax-4*r;
posy_min = 4*r; posy_max = ymax-4*r;
posx_vec = posx_min:4*r:posx_max;
posy_vec = posy_min:4*r:posy_max;
nstep = 3;
Nmax = length(posx_vec)*length(posy_vec); % max number of adhesion in the space
NmaxCount = ceil(Nmax/nstep);
expName = ['dia' num2str(2*r)];
rmsErrorL2 = zeros(nExp,NmaxCount);
rmsErrorL1 = zeros(nExp,NmaxCount);
EOA_L2 = zeros(nExp,NmaxCount);
EOA_L1 = zeros(nExp,NmaxCount);
lcurvePathL2 = cell(nExp,NmaxCount);
lcurvePathL1 = cell(nExp,NmaxCount);
fMapL2 = cell(nExp,NmaxCount);
fMapL1 = cell(nExp,NmaxCount);
oMapL2 = cell(nExp,NmaxCount);
oMapL1 = cell(nExp,NmaxCount);
bead_x = cell(nExp,1);
bead_y = cell(nExp,1);
Av = cell(nExp,1);
cropInfoL2 = cell(nExp,NmaxCount);
cropInfoL1 = cell(nExp,NmaxCount);
%% Simulations - simple enough simulation - only d
for epm=1:nExp
    ii=0;
    for n=1:nstep:Nmax % distance between adhesions
        ii=ii+1;
        if ii==1
            method = 'L2';
            dataPath=['/project/cellbiology/gdanuser/adhesion/Sangyoon/TFM simulations/AdhDensity/' expName '/exp' num2str(epm) 'f' num2str(f) 'n' num2str(n) method];
            [rmsErrorL2(epm,ii),EOA_L2(epm,ii),lcurvePathL2{epm,ii},fMapL2{epm,ii},cropInfoL2{epm,ii},oMapL2{epm,ii},bead_x{epm}, bead_y{epm}, Av{epm}] = ...
                testForceAdhDensity(d,f,r,n,method,dataPath);
            method = 'L1';
            oldDataPath = dataPath;
            dataPath=['/project/cellbiology/gdanuser/adhesion/Sangyoon/TFM simulations/AdhDensity/' expName '/exp' num2str(epm) 'f' num2str(f) 'n' num2str(n) method];
            [rmsErrorL1(epm,ii),EOA_L1(epm,ii),lcurvePathL1{epm,ii},fMapL1{epm,ii},cropInfoL1{epm,ii},oMapL1{epm,ii}] = ...
                testForceAdhDensity(d,f,r,n,method,dataPath,oldDataPath);
        else
            method = 'L2';
            dataPath=['/project/cellbiology/gdanuser/adhesion/Sangyoon/TFM simulations/AdhDensity/' expName '/exp' num2str(epm) 'f' num2str(f) 'n' num2str(n) method];
            [rmsErrorL2(epm,ii),EOA_L2(epm,ii),lcurvePathL2{epm,ii},fMapL2{epm,ii},cropInfoL2{epm,ii},oMapL2{epm,ii}] = ...
                testForceAdhDensity(d,f,r,n,method,dataPath,[],bead_x{epm}, bead_y{epm}, Av{epm});
            oldDataPath = dataPath;
            method = 'L1';
            dataPath=['/project/cellbiology/gdanuser/adhesion/Sangyoon/TFM simulations/AdhDensity/' expName '/exp' num2str(epm) 'f' num2str(f) 'n' num2str(n) method];
            [rmsErrorL1(epm,ii),EOA_L1(epm,ii),lcurvePathL1{epm,ii},fMapL1{epm,ii},cropInfoL1{epm,ii},oMapL1{epm,ii}] = ...
                testForceAdhDensity(d,f,r,n,method,dataPath,oldDataPath);
        end
    end
end

%% Plotting
