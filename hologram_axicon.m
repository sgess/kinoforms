% Number of points and size of mask
n_pts   = 6000;  % number of pixels
r_mask  = 20000; % radius of mask in microns
window  = 20;   % displayed region of mask/image 
z       = 5.0E5; % propagation distance in microns

% Ti:Saph parameters
lambda  = 0.800; % wavelength in microns
k_0     = 2*pi/lambda; % microns^-1

% Convergence angle
gamma = 0.50; % degrees
kp    = 2*pi*sind(gamma)/lambda; % microns^-1

% Mask coordinates
x  = linspace(-r_mask,r_mask,n_pts);
dx = x(2)-x(1);
xx = repmat(x,n_pts,1);
yy = rot90(xx,1);
r2 = xx.^2+yy.^2;
r  = sqrt(r2);
phi = atan2(yy,xx);
phi(round((n_pts/2+1)):n_pts,:) = phi(round((n_pts/2+1)):n_pts,:)+2*pi;

% Axis at image plane
fx = linspace(-lambda*z/(2*dx),lambda*z/(2*dx),n_pts);
dfx = fx(2)-fx(1);
Dfx = fx(end)-fx(1);

% Phase pattern
PSI = kp*r;

% Gaussian envelope (gets rid of FFT artifacts)
sig = 1;
gauss_env = exp(-r2/(2*(r_mask/sig)^2));

% Mask patterns
binary_mask = true(n_pts);
binary_mask(rem(PSI,2*pi) > pi) = false;
binary_mask = double(pi*binary_mask);

stair8_mask = zeros(n_pts);
stair8_mask(rem(PSI,2*pi) > 2*pi/8)  = 2*pi/8;
stair8_mask(rem(PSI,2*pi) > 4*pi/8)  = 4*pi/8;
stair8_mask(rem(PSI,2*pi) > 6*pi/8)  = 6*pi/8;
stair8_mask(rem(PSI,2*pi) > 8*pi/8)  = 8*pi/8;
stair8_mask(rem(PSI,2*pi) > 10*pi/8) = 10*pi/8;
stair8_mask(rem(PSI,2*pi) > 12*pi/8) = 12*pi/8;
stair8_mask(rem(PSI,2*pi) > 14*pi/8) = 14*pi/8;
stair8_mask(rem(PSI,2*pi) > 16*pi/8) = 16*pi/8;

% Field after mask
u_ideal  = gauss_env.*exp(-1i*PSI);
u_binary = gauss_env.*exp(-1i*binary_mask);
u_stair8 = gauss_env.*exp(-1i*stair8_mask);

% fresnel term
fresnel = exp(1i*k_0*r2/(2*z));

% fresnel diffraction
img_ideal = fft2(u_ideal.*fresnel);
img_ideal = fftshift(img_ideal);
i_ideal   = img_ideal.*conj(img_ideal);
i_ideal   = (k_0/z)^2*i_ideal;

img_binary = fft2(u_binary.*fresnel);
img_binary = fftshift(img_binary);
i_binary   = img_binary.*conj(img_binary);
i_binary   = (k_0/z)^2*i_binary;

img_stair8 = fft2(u_stair8.*fresnel);
img_stair8 = fftshift(img_stair8);
i_stair8   = img_stair8.*conj(img_stair8);
i_stair8   = (k_0/z)^2*i_stair8;

% Plots
mid = n_pts/2 + 1;
low = mid-window;
high = mid+window;
im_low = mid-window;
im_high = mid+window;

figure(1); subplot(2,2,1);
imagesc(fx(low:high),fx(low:high),PSI(low:high,low:high)); colormap gray; set(gca,'fontsize',14);
xlabel('X [\mum]','fontsize',16); ylabel('Y [\mum]','fontsize',16); title('Input Phase','fontsize',16);
subplot(2,2,2);
imagesc(fx(im_low:im_high),fx(im_low:im_high),i_ideal(im_low:im_high,im_low:im_high)); colormap gray; set(gca,'fontsize',14);
xlabel('X [\mum]','fontsize',16); ylabel('Y [\mum]','fontsize',16); title('Image Intensity','fontsize',16);
subplot(2,2,4);
plot(fx(im_low:im_high),i_ideal(im_low:im_high,mid)); axis tight; set(gca,'fontsize',14);
xlabel('X [\mum]','fontsize',16); ylabel('Intensity [arb]','fontsize',16); title('Lineout','fontsize',16);
set(gcf,'color','w');

figure(2); subplot(2,2,1);
imagesc(fx(low:high),fx(low:high),binary_mask(low:high,low:high)); colormap gray; set(gca,'fontsize',14);
xlabel('X [\mum]','fontsize',16); ylabel('Y [\mum]','fontsize',16); title('Input Phase','fontsize',16);
subplot(2,2,2);
imagesc(fx(im_low:im_high),fx(im_low:im_high),i_binary(im_low:im_high,im_low:im_high)); colormap gray; set(gca,'fontsize',14);
xlabel('X [\mum]','fontsize',16); ylabel('Y [\mum]','fontsize',16); title('Image Intensity','fontsize',16);
subplot(2,2,4);
plot(fx(im_low:im_high),i_binary(im_low:im_high,mid)); axis tight; set(gca,'fontsize',14);
xlabel('X [\mum]','fontsize',16); ylabel('Intensity [arb]','fontsize',16); title('Lineout','fontsize',16);
set(gcf,'color','w');

figure(3); subplot(2,2,1);
imagesc(fx(low:high),fx(low:high),stair8_mask(low:high,low:high)); colormap gray; set(gca,'fontsize',14);
xlabel('X [\mum]','fontsize',16); ylabel('Y [\mum]','fontsize',16); title('Input Phase','fontsize',16);
subplot(2,2,2);
imagesc(fx(im_low:im_high),fx(im_low:im_high),i_stair8(im_low:im_high,im_low:im_high)); colormap gray; set(gca,'fontsize',14);
xlabel('X [\mum]','fontsize',16); ylabel('Y [\mum]','fontsize',16); title('Image Intensity','fontsize',16);
subplot(2,2,4);
plot(fx(im_low:im_high),i_stair8(im_low:im_high,mid)); axis tight; set(gca,'fontsize',14);
xlabel('X [\mum]','fontsize',16); ylabel('Intensity [arb]','fontsize',16); title('Lineout','fontsize',16);
set(gcf,'color','w');