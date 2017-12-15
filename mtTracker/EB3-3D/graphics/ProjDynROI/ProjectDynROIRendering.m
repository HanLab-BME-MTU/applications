classdef ProjectDynROIRendering < ProjectDynROIProcess
%% This class provide a view from an orthogonal projection rawProjectDynROIProcess.
%% It reuses ProjectDynROIProcess as a backEnd, since this class already the necessary containers. 
    methods 

    function obj = ProjectDynROIRendering(rawProjectDynROIProcess,name)
        obj=obj@ProjectDynROIProcess(rawProjectDynROIProcess.getOwner(),'','nRenderedChannel',1);

        if(nargin>1)
            obj.buildAndSetOutFilePaths([rawProjectDynROIProcess.getOutputDir() filesep 'Rendering' filesep name],1);
            set(obj,'ref',rawProjectDynROIProcess.ref);
            set(obj,'nFrames',length(rawProjectDynROIProcess.nFrames));   
            [BX,BY,BZ]=rawProjectDynROIProcess.getBoundingBox();
            obj.setBoundingBox(BX,BY,BZ);
        end
        end

end
end
