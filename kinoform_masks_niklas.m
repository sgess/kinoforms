% Script to generate masks for kinoform phase plates

n_pts  = 12000;  % number of pixels
r_mask = 20000;  % radius of mask in microns
lambda = 0.800; % wavelength in microns
gma_05 = 0.5;   % angle of focused rays with respect to axis in degrees
gma_10 = 1.0;   % angle of focused rays with respect to axis in degrees

k_0    = 2*pi/lambda;              % microns^-1
kp_05  = 2*pi*sind(gma_05)/lambda; % microns^-1
kp_10  = 2*pi*sind(gma_10)/lambda; % microns^-1

% map x-y into r-phi
x   = linspace(-r_mask,r_mask,n_pts);
xx  = repmat(x,n_pts,1);
yy  = rot90(xx,1);
r   = sqrt(xx.^2+yy.^2);
phi = atan2(yy,xx);
% need to add 2pi to get correct phase in chosen coordinates
phi((n_pts/2+1):n_pts,:) = phi((n_pts/2+1):n_pts,:)+2*pi;

% Phase function PSI = k_perp*r + m*phi
PSI_m4kp05 = kp_05*r+4*phi;
PSI_m4kp10 = kp_10*r+4*phi;
PSI_m5kp05 = kp_05*r+5*phi;
PSI_m5kp10 = kp_10*r+5*phi;

% Generate masks and apply phase condition
mask_m4kp05 = false(n_pts);
mask_m4kp10 = false(n_pts);
mask_m5kp05 = false(n_pts);
mask_m5kp10 = false(n_pts);
mask_m4kp05(rem(PSI_m4kp05,2*pi) > pi) = true;
mask_m4kp10(rem(PSI_m4kp10,2*pi) > pi) = true; 
mask_m5kp05(rem(PSI_m5kp05,2*pi) > pi) = true; 
mask_m5kp10(rem(PSI_m5kp10,2*pi) > pi) = true; 


% Plot central 10% of mask (entire mask is way too big to plot)
x_ind  = abs(x)  < r_mask/10;

figure(1);
colormap gray;
imagesc(x(x_ind),x(x_ind),mask_m4kp05(x_ind,x_ind));
xlabel('microns','fontsize',14);
ylabel('microns','fontsize',14);
title('Kinoform Mask for J_4, \lambda = 800 nm, \gamma = 0.5^o','fontsize',14);

figure(2);
colormap gray;
imagesc(x(x_ind),x(x_ind),mask_m4kp10(x_ind,x_ind));
xlabel('microns','fontsize',14);
ylabel('microns','fontsize',14);
title('Kinoform Mask for J_4, \lambda = 800 nm, \gamma = 1^o','fontsize',14);

figure(3);
colormap gray;
imagesc(x(x_ind),x(x_ind),mask_m5kp05(x_ind,x_ind));
xlabel('microns','fontsize',14);
ylabel('microns','fontsize',14);
title('Kinoform Mask for J_5, \lambda = 800 nm, \gamma = 0.5^o','fontsize',14);

figure(4);
colormap gray;
imagesc(x(x_ind),x(x_ind),mask_m5kp10(x_ind,x_ind));
xlabel('microns','fontsize',14);
ylabel('microns','fontsize',14);
title('Kinoform Mask for J_5, \lambda = 800 nm, \gamma = 1^o','fontsize',14);