function likelihood = neg2LnLikelihoodXMEXcall(paramV,prob)

trajLength = size(prob.user.trajOut.observations,1);

prob.user.xOrder = length(paramV) - (prob.user.arOrder + prob.user.maOrder) - 1;

if prob.user.xOrder < 1
   prob.user.trajIn.observations = zeros(trajLength,2);   
end


prob.user.TRAJ = cat(3, prob.user.trajIn.observations , prob.user.trajOut.observations);
prob.user.TRAJ = permute( prob.user.TRAJ, [1 3 2]);
prob.user.wnVariance = 1;
prob = prob.user;

likelihood = neg2LnLikelihoodXMEX(paramV,prob);

