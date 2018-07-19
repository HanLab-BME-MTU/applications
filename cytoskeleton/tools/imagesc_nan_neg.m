function imagesc_nan_neg(inout_im,neg_value)

inout_im(isnan(inout_im)) = neg_value;

imagesc(inout_im);