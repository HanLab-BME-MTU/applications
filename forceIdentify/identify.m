%Set up monitor of the computation time.
startTime = cputime;

%Set up the model including the geometry, the mesh, the coefficient functions
% in the equation and the boundary condition.
setupModel;

%%%% Start the identification %%%%%%%%%

%Step 1: Construct the finite dimensional linear system to be solved for the 
% solution.
%Step 2 : Solve the linear system.
%Step 3 : Post processing and analysis of the results.
if strcmp(forceToIdentify,'bf') == 1
   calFwdOpBF;
   calRHVecBF;
   solveLSBF;
   postAssembleBF;
elseif strcmp(forceToIdentify,'tf') == 1
   calFwdOpTF;
   calRHVecTF;
   solveLSTF;
   postAssembleTF;
else
   error('Invalid value for ''forceToIdentify''.');
end

fprintf(1,'\nAll done : %d sec.\n', cputime-startTime);

