%Specify the method used to create the spline. FEMLAB provides four methods:
% 1. uniform 2. chordlength 3. centripetal 4. foley
splnMethod = 'chordlength';

%Load the geometry data of the cell.
%Get LEdgPx LEdgPy myoBFx myoBFy myoBBx myoBBy TEdgPx TEdgPy BEdgPx BEdgPy 
% ctrX ctrY smplPGx smplPGy bfDomPGx bfDomPGy recPGx recPGy;
load data/cellGeom;

c1 = geomspline([LEdgPx LEdgPy].','splinemethod',splnMethod,'closed','off');
c2 = geomspline([TEdgPx TEdgPy].','splinemethod',splnMethod,'closed','off');
c3 = geomspline([myoBFx myoBFy].','splinemethod',splnMethod,'closed','off');
c4 = geomspline([myoBBx myoBBy].','splinemethod',splnMethod,'closed','off');
rawGeom = geomcoerce('solid',{c1,c2,c3,c4});

%Rotate a little bit first so that the cell is horizontally placed.
angle = pi/10;
[LEdgPx,LEdgPy] = rotate2D(angle,LEdgPx-ctrX,-LEdgPy+ctrY);
[BEdgPx,BEdgPy] = rotate2D(angle,BEdgPx-ctrX,-BEdgPy+ctrY);
[TEdgPx,TEdgPy] = rotate2D(angle,TEdgPx-ctrX,-TEdgPy+ctrY);
[myoBFx,myoBFy] = rotate2D(angle,myoBFx-ctrX,-myoBFy+ctrY);
[myoBBx,myoBBy] = rotate2D(angle,myoBBx-ctrX,-myoBBy+ctrY);

c1 = geomspline([LEdgPx LEdgPy].','splinemethod',splnMethod,'closed','off');
c2 = geomspline([TEdgPx TEdgPy].','splinemethod',splnMethod,'closed','off');
c3 = geomspline([myoBFx myoBFy].','splinemethod',splnMethod,'closed','off');
c4 = geomspline([myoBBx myoBBy].','splinemethod',splnMethod,'closed','off');
geom = geomcoerce('solid',{c1,c2,c3,c4});

save data/simulGeom rawGeom geom ctrX ctrY angle;

