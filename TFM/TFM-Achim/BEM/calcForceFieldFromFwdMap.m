function [forceFieldNew]=calcForceFieldFromFwdMap(forceField,L)

forceFieldNew=forceField;

for i=1:length(forceField)
    M=forceField(i).par.M;
    u=forceField(i).par.u;
    forceMesh=forceField(i).par.forceMesh;
    [fx, fy, pos_x, pos_y]=calcSolFromFwdMapFastBEM(M,u,forceMesh,L);
    forceFieldNew(i).pos=horzcat(pos_x,pos_y); 
    forceFieldNew(i).vec=horzcat(fx,fy);
    % store all the paramters:
    forceFieldNew(i).par   = forceField(i).par;
    
    % Don't store M and insert the new regularization parameter:
    forceFieldNew(i).par.M =[];
    forceFieldNew(i).par.regParam=L;
end