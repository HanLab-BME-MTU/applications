function []=plotErrorForceSum(constrForceField)
errorForce=[];
for frame=1:length(constrForceField)
    if ~isempty(constrForceField{frame})
        try
            errorForce=vertcat(errorForce,constrForceField{frame}.errorSumForce.vec);
        catch
            % in case we give it the trackedNet structure:
            errorForce=vertcat(errorForce,constrForceField{frame}.stats.errorSumForce.vec);
        end
    end
end
plot(errorForce(:,1),errorForce(:,2),'-*');
axis equal
xlim([-max(abs(errorForce(:,1))) max(abs(errorForce(:,1)))])
ylim([-max(abs(errorForce(:,2))) max(abs(errorForce(:,2)))])
hold on
plot([0,0],[-max(abs(errorForce(:,2))) max(abs(errorForce(:,2)))],'--k')
plot([-max(abs(errorForce(:,1))) max(abs(errorForce(:,1)))],[0,0],'--k')