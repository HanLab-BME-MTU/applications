function imputed = myknnimpute(mat)
% myknnimpute Impute missing values in a matrix by using the most similar column
% vector (knnimpute.m). If all rows contain NaN, adjust to it. 
% 
% Jungsik Noh, 2016/01
% Updated, J Noh, 2017/08/28

% non-NaN proportion
cind = (sum(~isnan(mat)) ./ size(mat, 1)) > 0.95;

if sum(cind) > 0
    mat1 = mat(:, cind);

    % knnimpute.m
    immat1 = knnimpute(mat1);
    imputed = nan(size(mat));
    imputed(:, cind) = immat1;
else
    disp('myknnimpute failed and returned all NaNs.')
    imputed = nan(size(mat));
end

%imputed(:, ~cind) = mat(:, ~cind);

end
