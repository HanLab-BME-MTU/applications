%Post assemble of the identified boundary traction force and some calculation
% for debugging purpose. It is run after 'solveLSTF'.

fprintf(1,'Post assemble of the identified boundary traciton force : ');

localStartTime = cputime;

ll = 1;

%Assemble the identified boundary Traction force.
spTFx = cell(numEdges,1);
spTFy = cell(numEdges,1);
TFxR  = cell(numEdges,1);
TFyR  = cell(numEdges,1);

ampTFR = cell(numEdges,1);
phsTFR = cell(numEdges,1);
edgeFx = cell(numEdges,1);
edgeFy = cell(numEdges,1);

for k = 1:numEdges
   coefTF{k} = coef{k}(1:dimTF(k)).';
   spTFx{k}  = spmak(knotS{k},coefTF{k});

   coefTF{k} = coef{k}(dimTF(k)+1:2*dimTF(k)).';
   spTFy{k}  = spmak(knotS{k},coefTF{k});

   TFxR{k} = fnval(spTFx{k},edgeS{k});
   TFyR{k} = fnval(spTFy{k},edgeS{k});

   ampTFR{k} = sqrt(TFxR{k}.^2+TFyR{k}.^2);
   phsTFR{k} = angle(TFxR{k}+i*TFyR{k});
end

save([resultPath 'tfId'],'spTFx','spTFy');

%For debugging : Solve the elastic equation with the identified force and
% compare the computed displacement with the measured displacement data.

%Get the grid points that are inside the identification region.
[is,pe] = postinterp(fem,[gridPx gridPy].');
gridPx(pe) = [];
gridPy(pe) = [];

%Set the body force.
load([resultPath 'bfId']);
fn.BodyFx = [resultPath 'spMyoDFx'];
fn.BodyFy = [resultPath 'spMyoDFy'];
fp.BodyFx = {{'x' 'y'} {fs coefBF(1:end/2)}};
fp.BodyFy = {{'x' 'y'} {fs coefBF(end/2+1:end)}};

%Specify boundary condition.
for k = 1:numEdges
   fn.BndDispx{k} = 'bndDisp';
   fn.BndDispy{k} = 'bndDisp';
   fp.BndDispx{k} = {{'s'} {edgePPx{k}}};
   fp.BndDispy{k} = {{'s'} {edgePPy{k}}};

   fn.BndTracFx{k} = 'spBndTF';
   fn.BndTracFy{k} = 'spBndTF';
   fp.BndTracFx{k} = {{'s'} {spTFx{k}}};
   fp.BndTracFy{k} = {{'s'} {spTFy{k}}};
end

for k = 1:numEdges
   BCTypes{k} = 'Neumann';
   options    = elOptionsSet(options,'BCType',BCTypes);
   
   fem = elModelUpdate(fem,'options',options,'fn',fn,'fp',fp);
   fem = elasticSolve(fem,[]);

   %The displacements computed with the identified force. To be compared with 
   % 'dataU1' and 'dataU2'.
   [dataUC1 dataUC2] = postinterp(fem,'u1','u2',[dataPx dataPy].');
   [gridUC1 gridUC2] = postinterp(fem,'u1','u2',[gridPx gridPy].');

   [edgeUC1{k} edgeUC2{k}] = postinterp(fem,'u1','u2',edgeP{k});

   %The boundary traction force on the edge to be compared with 'TFxR' etc.
   [edgeFx{k} edgeFy{k}] = postinterp(fem,'g1','g2',edgeP{k});

   BCTypes{k} = 'Dirichlet';
end

fprintf(1,'%f sec.\n',cputime-localStartTime);
