function im_out = pad_boundary(im_in)

img_width = size(im_in,2);
img_height = size(im_in,1);
pad_xy = [img_width/2 img_height/2];

im_out = [ im_in im_in; im_in im_in;]*0;

im_out(pad_xy(2):pad_xy(2)+img_height-1, pad_xy(1):pad_xy(1)+img_width-1)=im_in;
