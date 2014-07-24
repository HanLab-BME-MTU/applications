function constrForceField=perfNetworkAnalysis(constrForceField,frame)

% determine the connectivity of each point along each interface:
constrForceField=createConctField(constrForceField,frame);

% now break the interface into smaller parts based on the connectivity
% of the points along each interface to determine interface-pieces
% that belong EXACTLY to two cells:
constrForceField=createTwoCellIntf(constrForceField,frame,0);

% Now create a force network:
constrForceField=createNetwork(constrForceField,frame);

% Propagate forces to find forces at the cell-cell interfaces. This
% only works through completely if there are no loops in the network:
constrForceField=propagateForces(constrForceField,frame);

%save([target_dir, filesep, 'cellCellForces.mat'], 'constrForceField', 'constrDisplField');
