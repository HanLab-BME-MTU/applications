function batchRun3

initialState = struct('mtLength0',0.44,'capSize0',5,'unitConc',10);

modelParam = struct('minLength',0.2,'maxLength',0.9','kTOnTFree',0.85,...
    'addAmpT',0.4,'addWidT',10,'addLenT',0.55,'kTOff',0.01,'kTOnD',0.45,...
    'kDOffFree',14,'addAmpD',7,'addWidD',10,'addLenD',0.55,'kHydrolysis',12);

runInfo = struct('maxNumSim',1,'totalTime',20000,'simTimeStep',0.01,...
    'timeEps',0.5,'expTimeStep',1,'aveInterval',0.6);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

hingeParam = struct('free',0,'diffConst',0);

hingeInit = struct('radius',0.1,'theta',pi/5,'phi',pi/3);

errFlag = analyzeHingeModel(2,modelParam,initialState,hingeParam,hingeInit,runInfo,'trial0000');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

hingeParam = struct('free',0,'diffConst',0.0001);

hingeInit = struct('radius',0.1,'theta',pi/5,'phi',pi/3);

errFlag = analyzeHingeModel(2,modelParam,initialState,hingeParam,hingeInit,runInfo,'trial0001');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

hingeParam = struct('free',0,'diffConst',0.0005);

hingeInit = struct('radius',0.1,'theta',pi/5,'phi',pi/3);

errFlag = analyzeHingeModel(2,modelParam,initialState,hingeParam,hingeInit,runInfo,'trial0005');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

hingeParam = struct('free',0,'diffConst',0.001);

hingeInit = struct('radius',0.1,'theta',pi/5,'phi',pi/3);

errFlag = analyzeHingeModel(2,modelParam,initialState,hingeParam,hingeInit,runInfo,'trial0010');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

hingeParam = struct('free',0,'diffConst',0.005);

hingeInit = struct('radius',0.1,'theta',pi/5,'phi',pi/3);

errFlag = analyzeHingeModel(2,modelParam,initialState,hingeParam,hingeInit,runInfo,'trial0050');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

hingeParam = struct('free',0,'diffConst',0.01);

hingeInit = struct('radius',0.1,'theta',pi/5,'phi',pi/3);

errFlag = analyzeHingeModel(2,modelParam,initialState,hingeParam,hingeInit,runInfo,'trial0100');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load trial0000;
dataStats11 = dataStats;
dataStats13 = dataStats3;

load trial0001;
dataStats21 = dataStats;
dataStats23 = dataStats3;

load trial0005;
dataStats31 = dataStats;
dataStats33 = dataStats3;

load trial0010;
dataStats41 = dataStats;
dataStats43 = dataStats3;

load trial0050;
dataStats51 = dataStats;
dataStats53 = dataStats3;

load trial0100;
dataStats61 = dataStats;
dataStats63 = dataStats3;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

gSpeed1 = [dataStats11.growthSpeed;dataStats21.growthSpeed;dataStats31.growthSpeed;dataStats41.growthSpeed;dataStats51.growthSpeed;dataStats61.growthSpeed];
gSpeed3 = [dataStats13.growthSpeed;dataStats23.growthSpeed;dataStats33.growthSpeed;dataStats43.growthSpeed;dataStats53.growthSpeed;dataStats63.growthSpeed];

sSpeed1 = [dataStats11.shrinkageSpeed;dataStats21.shrinkageSpeed;dataStats31.shrinkageSpeed;dataStats41.shrinkageSpeed;dataStats51.shrinkageSpeed;dataStats61.shrinkageSpeed];
sSpeed3 = [dataStats13.shrinkageSpeed;dataStats23.shrinkageSpeed;dataStats33.shrinkageSpeed;dataStats43.shrinkageSpeed;dataStats53.shrinkageSpeed;dataStats63.shrinkageSpeed];

catFreq1 = [dataStats11.catFreq;dataStats21.catFreq;dataStats31.catFreq;dataStats41.catFreq;dataStats51.catFreq;dataStats61.catFreq];
catFreq3 = [dataStats13.catFreq;dataStats23.catFreq;dataStats33.catFreq;dataStats43.catFreq;dataStats53.catFreq;dataStats63.catFreq];

resFreq1 = [dataStats11.resFreq;dataStats21.resFreq;dataStats31.resFreq;dataStats41.resFreq;dataStats51.resFreq;dataStats61.resFreq];
resFreq3 = [dataStats13.resFreq;dataStats23.resFreq;dataStats33.resFreq;dataStats43.resFreq;dataStats53.resFreq;dataStats63.resFreq];

save('tempData','gSpeed1','gSpeed3','sSpeed1','sSpeed3','catFreq1','catFreq3','resFreq1','resFreq3');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
