function [labels] = MRFEnergyMin_BeliefPropagation( Order1CliquePotential , Order2Cliques , Order2CliquePotential , varargin )
% 
% Input:
% 
% Order1_CliquePotential - <numSites x numLabels> matrix of real entries.
% Order2_Cliques - <numEdges x 2> full matrix of int entries.
% Order2_CliquePotential - <numLabels x numLabels x numEdges> full matrix of real entries.
% max_iters (vararg) - maximum iterations of the optimizer (generally 100-500 is sufficient)
% BP_TYPE (vararg) - 'LBP' for loopy belief propagation
%					 'TRBP' for tree-reweighted belief propagation
%
% Output:
% 
% flow - maximum flow value
% labels - a vector of size Nx1 containing the label of each node 
%          respectively
% 

p = inputParser;

p.addParamValue( 'max_iters' , 100 , @(x)isnumeric(x) && numel(x) == 1 && x > 0 );
p.addParamValue( 'BP_TYPE' , 'TRBP' , @(x)ischar(x) && ( strcmp( x , 'LBP' ) || strcmp( x , 'TRBP' ) ) );
p.parse( varargin{:} );

%% Optimize MRF using Loopy Belief Propagation

start_time = datevec( datestr( now ) );

fprintf( 1 , '\n\nPre-processing ... ' );

nNodes = size( Order1CliquePotential , 1 );
nStates = size( Order1CliquePotential , 2 );

nodePot = exp( -1 * Order1CliquePotential );
edgePot = exp( -1 * Order2CliquePotential );

edgeStruct.edgeEnds = Order2Cliques;
[ edgeStruct.V , edgeStruct.E ] = UGM_makeEdgeVE( Order2Cliques , nNodes );
edgeStruct.nNodes = nNodes;
edgeStruct.nEdges = size( Order2Cliques , 1 );
edgeStruct.nStates = nStates * ones( nNodes , 1 );
edgeStruct.useMex = 1;
edgeStruct.maxIter = p.Results.max_iters;

fprintf( 1 , '\n\nRunning Belief Propagation ... ' );

if strcmp( p.Results.BP_TYPE , 'LBP' )

	labels = UGM_Decode_LBP(nodePot,edgePot,edgeStruct);
	
else
	
	labels = UGM_Decode_TRBP(nodePot,edgePot,edgeStruct);
	
end

end_time = datevec( datestr( now ) );

fprintf( 1 , '\n\nMRF Optimization took %f seconds ... \n\n' ,  etime( end_time , start_time ) );

clear nodePot edgePot edgeStruct;

end