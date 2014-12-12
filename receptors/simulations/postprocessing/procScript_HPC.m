%   HPC script to generate cluster statistics for all simulations. This
%   script is used for short time step runs (a millisecond or smaller) 
%   which makes calculating cluster statistics on the desktop impractical
%   or impossible because of the amount of memory required.
%   
%   The root directory, recept2clustAssign and simParam mat files need to
%   be specified below. A shell script excutes this script and the cluster
%   statistics will saved to a mat file.
%
%   Robel Yirdaw, 04/03/14
%

try
    cd('/project/biophysics/jaqaman_lab/interKinetics/ryirdaw/2014/08/082614/dR0p25_aP1_dT0p001_all');
    
    load('recept2clustAssignAll_dR0p25_aP1_dT0p001_all.mat');
    load('valuesFor_dR0p25_aP1_dT0p001_all.mat', 'simParamAll');
    
    clusterStat = calcClusterStatsFromR2C(recept2clustAssignAll,simParamAll{1}.observeSideLen,simParamAll{1}.probDim);

    fprintf('Done creating clusterStat.  \nNow saving...');
    save('clusterStat_dR0p25_aP1_dT0p001_all.mat','clusterStat','-v7.3');
    
    fprintf('done.');
catch
    
end
