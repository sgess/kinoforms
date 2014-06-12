%win_fact = 1;

n_pts  = 750;  % number of pixels
r_mask = 20000;  % radius of mask in microns
lambda = 0.800;  % wavelength in microns
gma_300 = 0.312;   % angle of focused rays with respect to axis in degrees
z = 5E5;

k_0    = 2*pi/lambda;              % microns^-1
kp_300  = 2*pi*sind(gma_300)/lambda; % microns^-1

x   = linspace(-r_mask,r_mask,n_pts);
[xx,yy]  = meshgrid(x,x);
r2  = xx.^2+yy.^2;
r = sqrt(r2);
%phi = atan2(yy,xx);
% need to add 2pi to get correct phase in chosen coordinates
%phi(round((n_pts/2+1)):n_pts,:) = phi(round((n_pts/2+1)):n_pts,:)+2*pi;

sig = 4;
gauss_env = exp(-r2/(2*(r_mask/sig)^2));
%mask_env = r < r_mask/sig;
mask_env = ones(n_pts);
smooth_env = gauss_env.*mask_env;

PSI_axi = kp_300*r;
%mask_axi = false(n_pts);
%mask_axi(rem(PSI_axi,2*pi) > pi) = true;

%u_axi = smooth_env.*exp(-1i*mask_axi);
u_axi = smooth_env.*exp(-1i*PSI_axi);
fresnel = exp(1i*k_0*r2/(2*z));

%i_axi = u_axi.*conj(u_axi);
%figure(1); imagesc(i_axi); colorbar;

%u_img = fft2(u_axi);
u_img = fft2(u_axi.*fresnel);
u_img = fftshift(u_img);
u_img = exp(1i*k_0*z)*u_img/(1i*lambda*z);
i_img = u_img.*conj(u_img);
%figure(2); imagesc(log(i_img+1)); colorbar; colormap gray;
figure(2); imagesc(i_img); colorbar; colormap gray;
