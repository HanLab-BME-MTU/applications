function [randTracksCell]=randomizeTracksMC(MD,randomDist,varargin)
ip = inputParser;
ip.CaseSensitive = false;
ip.KeepUnmatched = true;
ip.addRequired('MD',@(MD) isa(MD,'MovieData'));
ip.addRequired('randomDist');
ip.addOptional('dynManifoldsCell',[]);
ip.addParameter('randomType','manifoldAntispace');
ip.addParameter('tracks',[]);
ip.addParameter('mappingDistance',10);
ip.addParameter('simuNumber',100);
ip.addParameter('process',[]);
ip.parse(MD,randomDist,varargin{:});
p=ip.Results;


randTracksCell=cell(1,p.simuNumber)

for sIdx=1:p.simuNumber
    randTracksCell{sIdx}=randomizeTracksMC(MD,randomDist,varargin{:});
end
process=p.process;
if(~isempty(process))
    %%
    procFolder=[process.getOwner().outputDirectory_  filesep 'Kin' filesep 'randomized' filesep];
    mkdirRobust(procFolder);
    save([procFolder 'MC-' num2str(p.simuNumber) p.randomType '-randKinTracks.mat'],'randTracksCell');
    process.setOutFilePaths({[procFolder p.randomType '-randKinTracks.mat']})
    pa = process.getParameters();
    pa = ip.Results;
    process.setParameters(pa);
    process.setDateTime();
end;

end
