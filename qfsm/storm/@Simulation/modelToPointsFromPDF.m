function modelToPointsFromPDF(obj,model,sigmaNoise,pdf)
if ~isempty(model)
    obj.data.points = [];
    for m=1:numel(model)
        obj.filamentFromPDF(model{m},sigmaNoise,pdf);
    end
    obj.data.rawPoints = obj.data.points;
    disp('Simulation: Points generated from the model!');
else
    disp('Simulation: No model!');
end
end