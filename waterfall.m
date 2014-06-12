% input parameters 
n_pts   = 3000;              % number of pixels
r_mask  = 20000;             % radius of mask in microns
lambda  = 0.800;             % wavelength in microns
gma_300 = 0.345;             % angle of focused rays [degrees]
k_0     = 2*pi/lambda;       % microns^-1
kp_300  = k_0*sind(gma_300); % microns^-1
z_start = 5E5;               % plasma start
z_end   = 4E6;               % plasma end
n_z     = 351;               % plasma points
laser_r = 20000;             % laser radius, microns
laser_E = 0.500;             % laser energy, J
laser_t = 100e-15;           % laser pulse width, s

laser_I0 = laser_E/(laser_t*pi*(laser_r/1e4)^2);


% initialize grid, n_pts x n_pts
x = linspace(-r_mask,r_mask,n_pts);
[xx,yy]  = meshgrid(x,x);
r2 = xx.^2+yy.^2;
r = sqrt(r2);

%phi = atan2(yy,xx);
% need to add 2pi to get correct phase in chosen coordinates
%phi(round((n_pts/2+1)):n_pts,:) = phi(round((n_pts/2+1)):n_pts,:)+2*pi;

% intensity envelope
sig = 0.95;
%gauss_env = exp(-r2/(2*(r_mask/sig)^2));
gauss_env = ones(n_pts);
mask_env = r < sig*r_mask;
%mask_env = ones(n_pts);
smooth_env = gauss_env.*mask_env;

% initialize phase term
PSI_axi = kp_300*r;
%u_axi = laser_I0*smooth_env.*exp(-1i*PSI_axi);

% initialize mask phase
mask_axi = zeros(n_pts);
mask_axi(rem(PSI_axi,2*pi) > pi) = pi;
u_axi = laser_I0*smooth_env.*exp(-1i*mask_axi);

% loop over z
z = linspace(z_start,z_end,n_z);
water = zeros(n_pts,n_z);

for i = 1:n_z
    
    display(i);
    
    fresnel = exp(1i*k_0*r2/(2*z(i)));

    u_img = fft2(u_axi.*fresnel)/n_pts;
    u_img = fftshift(u_img);
    u_img = exp(1i*k_0*z(i))*u_img/(1i*lambda*z(i));
    i_img = u_img.*conj(u_img);
    
    water(:,i) = i_img(:,n_pts/2+1);

end

figure(1);
pcolor(z/1e6,x/1e3,water); shading flat; colormap hot; h = colorbar; caxis([0 3e12]);
xlabel('Z [m]','fontsize',16); ylabel('X [mm]','fontsize',16); 
ylabel(h,'I [W/cm^2]','fontsize',16); title('Laser Intensity','fontsize',16);
set(gca,'fontsize',16); set(gcf,'color','w'); saveas(gcf,'kino_2d_hr2.eps','epsc');

figure(2);
plot(z/1e6,water(n_pts/2+1,:),'linewidth',2);
xlabel('Z [m]','fontsize',16); ylabel('I [W/cm^2]','fontsize',16);
title('On-Axis Laser Intensity','fontsize',16);
set(gca,'fontsize',16); set(gcf,'color','w'); saveas(gcf,'kino_1d_hr2.eps','epsc');
