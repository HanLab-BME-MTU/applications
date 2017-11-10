function renderGEFDiffSimulation(MD,XY,time,diffCoef,diffTime)%,,diffCoeff,varargin)

activationMap=cell(1,MD.nFrames_);
activationMap=cellfun(@(x) zeros(MD.imSize_),activationMap,'unif',false);

for i=1:length(XY)
    for t=time(i):min(MD.nFrames_,(time(i)+diffTime))
        idx=t-time(i)+1;
        activationMap{t}=MidpointCircle(activationMap{t},diffCoef*idx+1,XY(i,2), XY(i,1),diffTime-idx);
    end
end
%%
mkdir([MD.outputDirectory_ filesep 'MT-GEF' filesep  'activationSimulation' filesep 'activationMap']);
parfor i=1:(MD.nFrames_)
   imwrite(uint8(activationMap{i}),[MD.outputDirectory_ filesep 'MT-GEF' filesep 'activationSimulation' filesep 'activationMap/map' num2str(i,'%04d') '.tif']);
end