% Script to generate masks for kinoform phase plates

n_pts  = 12000;  % number of pixels
r_mask = 20000;  % radius of mask in microns
lambda = 0.800;  % wavelength in microns
gma_300 = 0.312;   % angle of focused rays with respect to axis in degrees
gma_400 = 0.273;   % angle of focused rays with respect to axis in degrees
gma_500 = 0.250;   % angle of focused rays with respect to axis in degrees
gma_600 = 0.235;   % angle of focused rays with respect to axis in degrees

k_0    = 2*pi/lambda;              % microns^-1
kp_300  = 2*pi*sind(gma_300)/lambda; % microns^-1
kp_400  = 2*pi*sind(gma_400)/lambda; % microns^-1
kp_500  = 2*pi*sind(gma_500)/lambda; % microns^-1
kp_600  = 2*pi*sind(gma_600)/lambda; % microns^-1

% map x-y into r-phi
x   = linspace(-r_mask,r_mask,n_pts);
xx  = repmat(x,n_pts,1);
yy  = rot90(xx,1);
r   = sqrt(xx.^2+yy.^2);
phi = atan2(yy,xx);
% need to add 2pi to get correct phase in chosen coordinates
phi((n_pts/2+1):n_pts,:) = phi((n_pts/2+1):n_pts,:)+2*pi;

% Phase function PSI = k_perp*r + m*phi
PSI_300 = kp_300*r+5*phi;
PSI_400 = kp_400*r+6*phi;
PSI_500 = kp_500*r+7*phi;
PSI_600 = kp_600*r+8*phi;

% Generate masks and apply phase condition
mask_300 = false(n_pts);
mask_400 = false(n_pts);
mask_500 = false(n_pts);
mask_600 = false(n_pts);
mask_300(rem(PSI_300,2*pi) > pi) = true;
mask_400(rem(PSI_400,2*pi) > pi) = true;
mask_500(rem(PSI_500,2*pi) > pi) = true;
mask_600(rem(PSI_600,2*pi) > pi) = true;

% Plot central 10% of mask (entire mask is way too big to plot)
x_ind  = abs(x)  < r_mask/10;

figure(1);
colormap gray;
imagesc(x(x_ind),x(x_ind),mask_300(x_ind,x_ind));
xlabel('microns','fontsize',14);
ylabel('microns','fontsize',14);
title('Kinoform Mask for 300 \mum channel, J_5, \gamma = 0.312^o','fontsize',14);

figure(2);
colormap gray;
imagesc(x(x_ind),x(x_ind),mask_400(x_ind,x_ind));
xlabel('microns','fontsize',14);
ylabel('microns','fontsize',14);
title('Kinoform Mask for 400 \mum channel, J_6, \gamma = 0.273^o','fontsize',14);

figure(3);
colormap gray;
imagesc(x(x_ind),x(x_ind),mask_500(x_ind,x_ind));
xlabel('microns','fontsize',14);
ylabel('microns','fontsize',14);
title('Kinoform Mask for 500 \mum channel, J_7, \gamma = 0.250^o','fontsize',14);

figure(4);
colormap gray;
imagesc(x(x_ind),x(x_ind),mask_600(x_ind,x_ind));
xlabel('microns','fontsize',14);
ylabel('microns','fontsize',14);
title('Kinoform Mask for 600 \mum channel, J_8, \gamma = 0.235^o','fontsize',14);
