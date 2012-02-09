function modelToPoints(obj,model,sigmaNoise,density)
if ~isempty(model)
    obj.data.points = [];
    for m=1:numel(model)
        obj.bezierWithDensity(model{m},density);
        % ===
        % nSamples = ceil(lengthBezier(model{m})*density);
        % sourceSpacing = 5;
        % fractionOfDamagedSources = 0.5;
        % obj.filament(model{m},nSamples,sigmaNoise,sourceSpacing,fractionOfDamagedSources)
        % ===
    end
    obj.addGaussianNoise(sigmaNoise);
    obj.data.rawPoints = obj.data.points;
    disp('Simulation: Points generated from the model!');
else
    disp('Simulation: No model!');
end
end