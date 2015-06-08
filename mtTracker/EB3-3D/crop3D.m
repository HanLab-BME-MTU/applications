function crop3D(MD,ROI)    
MD.addProcess(ExternalProcess(MD))
MD.getProcess(1).setInFilePaths({})
MD.getProcess(1).setOutFilePaths({[MD.outputDirectory_ filesep 'cropped']});
p = MD.getProcess(1).getParameters();
p.metadata1 = 'crop3D';
p.metadata2 = ROI;
MD.getProcess(1).setParameters(p);

for cIdx=1:length(MD.channels_)
    mkdir([MD.getProcess(1).outFilePaths_{1} filesep 'ch' num2str(cIdx-1)]);
    for t=1:MD.nFrames_
        vol=MD.getChannel(cIdx).loadStack(t);
        stackWrite(vol(ROI(1):ROI(4),ROI(2):ROI(5),ROI(3):ROI(6)) ,[MD.getProcess(1).outFilePaths_{1} filesep 'ch' num2str(cIdx-1) filesep 'time-' num2str(t,'%04d') '.tif']);
    end 
end 
        
    
