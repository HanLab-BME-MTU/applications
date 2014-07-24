function [idBestMatch]=findBasisClassInTbl(basisClassTbl,basisClassIn,xrange,yrange,meshPtsFwdSol)


%**************************************************************************
% 1) Try to find the basis function in the tablebase:                     *
%**************************************************************************
[foundClass]=findBasisClass(basisClassTbl,basisClassIn.neighPos);

idUsable=[];
smplDnst=[];
idBestMatch=[];
% Check in all classes if the range and the meshPtsFwdSol fits:
for idClass=foundClass'
    % read out the important values:
    meshPtsFwdSolTbl=basisClassTbl(idClass).uSol.meshPtsFwdSol;
    xrangeTbl=basisClassTbl(idClass).uSol.xrange;
    yrangeTbl=basisClassTbl(idClass).uSol.yrange;
    
    check1 = (xrangeTbl(1)<=xrange(1) && xrangeTbl(2)>=xrange(2));
    check2 = (yrangeTbl(1)<=yrange(1) && yrangeTbl(2)>=yrange(2));
    check3 = (meshPtsFwdSolTbl>=meshPtsFwdSol);
    
    if check1 && check2 && check3
        % then we can use this tabled solution for this basis class
        idUsable=horzcat(idUsable,idClass);
        
        % store the sample density to select the best one:
        smplDnst=horzcat(smplDnst,meshPtsFwdSolTbl/((xrangeTbl(2)-xrangeTbl(1))*(yrangeTbl(2)-yrangeTbl(1))));
    end
end
if ~isempty(idUsable)
    % The best stored solution is the one with the highes smplDensity:
    [~,posBestMatch]=max(smplDnst);
    idBestMatch=idUsable(posBestMatch);
end