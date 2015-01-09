function [fn,fp,BCTypes,bndInd] = initializeCMBoundaryCondition(numEdges,borderE,borderSeg,bndInd)

      %Specify the boundary traction force.
      fn.BndDispx  = cell(1,numEdges);
      fn.BndDispy  = cell(1,numEdges);
      fp.BndDispx  = cell(1,numEdges);
      fp.BndDispy  = cell(1,numEdges);
      fn.BndTracFx = cell(1,numEdges);
      fn.BndTracFy = cell(1,numEdges);
      fp.BndTracFx = cell(1,numEdges);
      fp.BndTracFy = cell(1,numEdges);
      if ~isempty(borderE)
         fp.BndTracFy = cell(1,numEdges+1);
      else
         fp.BndTracFy = cell(1,numEdges);
      end
          
      for k = 1:numEdges
         %Specify zero boundary displacement initially.
         fn.BndDispx{k} = 0;
         fn.BndDispy{k} = 0;

         %Set up link to boundary traction force function.
         fn.BndTracFx{k} = 'spBndTF';
         fn.BndTracFy{k} = 'spBndTF';
         fp.BndTracFx{k} = {{'s'} {0}};
         fp.BndTracFy{k} = {{'s'} {0}};

         %Specify the type of boundary conditions.
         BCTypes{k} = 'Dirichlet';
      end

      if ~isempty(borderE)
         bndInd{numEdges+1} =  borderSeg;

         fn.BndDispx{numEdges+1}  = [];
         fn.BndDispy{numEdges+1}  = [];
         fp.BndDispx{numEdges+1}  = [];
         fp.BndDispy{numEdges+1}  = [];

         fn.BndTracFx{numEdges+1} = 0;
         fn.BndTracFy{numEdges+1} = 0;
         fp.BndTracFx{numEdges+1} = [];
         fp.BndTracFy{numEdges+1} = [];

         BCTypes{numEdges+1} = 'Neumann';
      end

      %Specify the body force.
      fn.BodyFx = 'femBodyF';
      fn.BodyFy = 'femBodyF';
      fp.BodyFx = {{'x' 'y'} {[] 0}};
      fp.BodyFy = {{'x' 'y'} {[] 0}};
      fn.YModul = 1; %
      fn.PRatio = 0.01; %
      fp.YModul = []; %
      fp.PRatio = []; % These need to be updated later.
%       fn.YModul = 'elImgParFun';
%       fp.YModul = {{'x' 'y'} {actImg 1.5 2 0.1 'Gaussian' 5}};
      

