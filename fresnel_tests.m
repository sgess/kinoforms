n_pts  = 12000;  % number of pixels
r_mask = 20000;  % radius of mask in microns
lambda = 0.800;  % wavelength in microns
gma_300 = 0.312;   % angle of focused rays with respect to axis in degrees

k_0    = 2*pi/lambda;              % microns^-1
kp_300  = 2*pi*sind(gma_300)/lambda; % microns^-1

x   = linspace(-r_mask,r_mask,n_pts);
dx = x(2)-x(1);
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

E0=1;
z = lambda*1E7;
U0=E0*exp(1i*PSI_300);
%U0=E0*exp(1i*mask_300);
ker=exp((1i*k_0/(2*z))*r2);

%%
x_ind  = abs(x)  < r_mask/10;
imagesc(x(x_ind),x(x_ind),PSI_300(x_ind,x_ind));
xlabel('X [\mum]','fontsize',16);
ylabel('Y [\mum]','fontsize',16);
title('\Phi = m\phi + k_{\perp}r','fontsize',16);
set(gca,'fontsize',16);

%%

perf = besselj(5,kp_300*r);
perf2 = perf.*perf;

imagesc(x(x_ind),x(x_ind),perf2(x_ind,x_ind));
xlabel('X [\mum]','fontsize',16);
ylabel('Y [\mum]','fontsize',16);
title('I \propto J_5^2(k_{\perp}r)','fontsize',16);
set(gca,'fontsize',16);

%%

imagesc(x(x_ind),x(x_ind),mask_300(x_ind,x_ind));
colormap gray;
xlabel('X [\mum]','fontsize',16);
ylabel('Y [\mum]','fontsize',16);
title('Binary phase grating','fontsize',16);
set(gca,'fontsize',16);



%%

U_twid = U0.*ker;
U = fft2(U_twid);
U2 = fftshift(U);
I = U2.*conj(U2);

%%
Unew = fresnel_advance(U0,dx,dx,z,lambda);
Unew2 = Unew.*conj(Unew);