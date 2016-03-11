function  GCAVisualsMakeOverlaysVeilStemWindows(windowC,normalsC,protC,edgeC)
%GCAVisualsMakeOverlaysWindows: 

%% INPUT 


%% OUTPUT
 if ~isempty(windowC{1}) ; % not sure why but it appears hunter leaves this empty sometimes ? 
      
 windowFirst = windowC{1}(1) ; 
 else 
     windowFirst = windowC{2}(1); 
     
 end 
     
  


            gcaPlotWindows(windowFirst,{'y','FaceAlpha',1},0,'bandMax',1); 
            gcaPlotWindows(windowC,{'g','FaceAlpha',0},5,'bandMax',1);   
            %Insure all normals for this frame have unit length - for some reason
            %Sam's software doesn't return them with unit length...
            magNorm = sqrt(dot(normalsC',normalsC')');
            normalsC = normalsC./ repmat(magNorm,1,2);
            
            %Get the normal component of the protrusion vectors for this frame
            protNorm = dot(protC',normalsC')';
            %And the magnitude of the protrusion vectors
            protMag = sqrt(dot(protC',protC'))';
            
            idxRetract = protNorm<0 & protNorm < - .462;
            
            idxProt = protNorm >0 & protNorm > .462;
            idxZero = ~(idxRetract | idxProt);
           
             quiver(edgeC(idxRetract,1),edgeC(idxRetract,2),protC(idxRetract,1),protC(idxRetract,2),0,'b');
             quiver(edgeC(idxZero,1), edgeC(idxZero,2), protC(idxZero,1),protC(idxZero,2),0,'g');
             quiver(edgeC(idxProt,1),edgeC(idxProt,2),protC(idxProt,1),protC(idxProt,2),0,'r');

end

