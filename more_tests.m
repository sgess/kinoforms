n_pts  = 12000;  % number of pixels
r_mask = 20000;  % radius of mask in microns
lambda = 0.800;  % wavelength in microns
gma_300 = 0.312;   % angle of focused rays with respect to axis in degrees

k_0    = 2*pi/lambda;              % microns^-1
kp_300  = 2*pi*sind(gma_300)/lambda; % microns^-1

x   = linspace(-r_mask,r_mask,n_pts);
xx  = repmat(x,n_pts,1);
yy  = rot90(xx,1);
r   = sqrt(xx.^2+yy.^2);
r2  = xx.^2+yy.^2;
phi = atan2(yy,xx);
% need to add 2pi to get correct phase in chosen coordinates
phi((n_pts/2+1):n_pts,:) = phi((n_pts/2+1):n_pts,:)+2*pi;

% Phase function PSI = k_perp*r + m*phi
PSI_300 = kp_300*r+5*phi;
mask_300 = false(n_pts);
mask_300(rem(PSI_300,2*pi) > pi) = true;
%u_300 = exp(-r2/5000).*exp(-1i*PSI_300);

circ_mask = false(n_pts);
circ_mask(r < 19999) = true;
u_300 = circ_mask.*exp(-1i*PSI_300);

% Phase function PSI = k_perp*r
PSI_axi = kp_300*r;
mask_axi = false(n_pts);
mask_axi(rem(PSI_axi,2*pi) > pi) = true;
u_axi = circ_mask.*exp(-1i*PSI_axi);

%%

ufft = fftshift(fft2(u_axi));
plot_u = ufft.*conj(ufft);

imagesc(log(plot_u));
colormap gray;
colorbar;
