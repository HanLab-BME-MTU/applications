function wm=wstarck222(x,k,irecon)

% Calculates the multiscale product from the reconstructed image 

[ly,lx]=size(x);		% Get the original image size


%------------------------------
%   SETUP ARRAYS 
%------------------------------
yl=zeros(ly,lx,k);				

tx=x;				% Copy the original.
yl(:,:,1)=x;

wm=ones(ly,lx);

     
  for ind=1:k			% For every scale...


	
	%%% Now let's call wtlo2 (low pass transform) first over the rows %%%
	%%% of the present image			   %%%

	
		tw(1:ly,1:lx)=wtlo2(tx(1:ly,1:lx),ind);
	

	%%% And then over the columns  %%%

	
		tw(1:ly,1:lx)=wtlo2(tw(1:ly,1:lx)',ind)';
	

	yl(:,:,ind+1)=tw;	% This is the ind_th approximation of x
			
	
    yh(:,:)=yl(:,:,ind)-yl(:,:,ind+1); % This is the ind_th detail or residue
    wm=abs(wm.*yh);
    %yh=tx-tw; % This is the ind_th detail or residue (or high pass component)
    tx=tw;   % This is the ind_th approximation of x (or low pass component)   
    st=std(yh(:)); % Estimate the standard deviation of the noise in the detail at each scale.
  end				% end of all scales.
   
