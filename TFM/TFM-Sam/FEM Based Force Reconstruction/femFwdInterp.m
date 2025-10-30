function vOut=femFwdInterp(x,y,allPts,dtIn,vIn,method)
    F=TriScatteredInterp(dtIn,vIn,method);
    vOut = F(x,y);
    checkVec = isnan(vOut);
    vOut(checkVec) = 0;
end