n_pts  = 12000;  % number of pixels
r_mask = 20000;  % radius of mask in microns
lambda = 0.800;  % wavelength in microns
gma_300 = 0.312;   % angle of focused rays with respect to axis in degrees

k_0    = 2*pi/lambda;              % microns^-1
kp_300  = 2*pi*sind(gma_300)/lambda; % microns^-1

r   = linspace(0,r_mask,n_pts/2);
z = lambda*5E6;
k   = k_0*r/z;
dr = r(2)-r(1);

%%

h = exp(1i*(kp_300*r+k_0*r.^2/z));

[H,I]=hat(h,r,k,5);
U = H.*conj(H);
ind = 1:600;
figure;
plot(r(ind),U(ind));

%%

h = exp(1i*(kp_300*r));

[H,I]=hat(h,r,k,5);
U = H.*conj(H);
ind = 1:600;
figure;
plot(r(ind),U(ind));