function []=plotErrorForceSum(constrForceField)
errorForce=[];
errorForce_perpix=[];
for frame=1:length(constrForceField)
    if ~isempty(constrForceField{frame})
        try
            errorForce        = vertcat(errorForce,constrForceField{frame}.errorSumForce.vec);
            area              = sum(constrForceField{frame}.segmRes.maskDilated(:));
            errorForce_perpix = vertcat(errorForce_perpix,constrForceField{frame}.errorSumForce.vec/area);
        catch
            % in case we give it the trackedNet structure:
            errorForce=vertcat(errorForce,constrForceField{frame}.stats.errorSumForce.vec);
            try
                errorForce_perpix=vertcat(errorForce_perpix,constrForceField{frame}.stats.errorSumForce.vec_perpix);
            end
        end
    end
end
figure()
plot(errorForce(:,1),errorForce(:,2),'-*b');
axis equal
xlim([-max(abs(errorForce(:,1))) max(abs(errorForce(:,1)))])
ylim([-max(abs(errorForce(:,2))) max(abs(errorForce(:,2)))])
hold on
plot([0,0],[-max(abs(errorForce(:,2))) max(abs(errorForce(:,2)))],'--k')
plot([-max(abs(errorForce(:,1))) max(abs(errorForce(:,1)))],[0,0],'--k')
xlabel('x-component')
ylabel('y-component')
title('Error force')

try
    errorForce_perpix(2,:);
    figure()
    plot(errorForce_perpix(:,1),errorForce_perpix(:,2),'-*r');
    axis equal
    xlim([-max(abs(errorForce_perpix(:,1))) max(abs(errorForce_perpix(:,1)))])
    ylim([-max(abs(errorForce_perpix(:,2))) max(abs(errorForce_perpix(:,2)))])
    hold on
    plot([0,0],[-max(abs(errorForce_perpix(:,2))) max(abs(errorForce_perpix(:,2)))],'--k')
    plot([-max(abs(errorForce_perpix(:,1))) max(abs(errorForce_perpix(:,1)))],[0,0],'--k')
    xlabel('x-component')
    ylabel('y-component')
    title('Error force per pixel')
end