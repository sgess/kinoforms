% Script to generate masks for kinoform phase plates

n_pts  = 12000;  % number of pixels
r_mask = 20000;  % radius of mask in microns
lambda = 0.800; % wavelength in microns
gma_300 = 0.312;   % angle of focused rays with respect to axis in degrees
gma_05  = 0.5;

k_0    = 2*pi/lambda;              % microns^-1
kp_300  = 2*pi*sind(gma_300)/lambda; % microns^-1
kp_05  = 2*pi*sind(gma_05)/lambda; % microns^-1

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
PSI_05 = kp_05*r+5*phi;
PSI_combo = PSI_05+PSI_300;

% Generate masks and apply phase condition
mask_300 = false(n_pts);
mask_05 = false(n_pts);
mask_combo = false(n_pts);

mask_300(rem(PSI_300,2*pi) > pi) = true;
mask_05(rem(PSI_05,2*pi) > pi) = true;
mask_combo(rem(PSI_combo,2*pi) > pi) = true;


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
imagesc(x(x_ind),x(x_ind),mask_05(x_ind,x_ind));
xlabel('microns','fontsize',14);
ylabel('microns','fontsize',14);
title('Kinoform Mask for 187 \mum channel, J_5, \gamma = 0.5^o','fontsize',14);

figure(3);
colormap gray;
imagesc(x(x_ind),x(x_ind),mask_combo(x_ind,x_ind));
xlabel('microns','fontsize',14);
ylabel('microns','fontsize',14);
title('Kinoform Mask for J_5 combination','fontsize',14);

figure(4);
colormap gray;
imagesc(x(x_ind),x(x_ind),PSI_combo(x_ind,x_ind));
xlabel('microns','fontsize',14);
ylabel('microns','fontsize',14);
title('Phase pattern for J_5 combination','fontsize',14);

% %
yi = linspace(0,300,301);
Intensity = zeros(size(yi));
z = 10000;

for i = 1:length(yi);
    
    sqr1 = xx.^2+(yy-yi(i)).^2;
    gamma = atan(sqrt(sqr1)/z);
    s = sqrt(z^2+sqr1);
    PHI=exp(1i*(PSI_combo+k_0*s));
    %PHI=exp(1i*(pi*double(mask_m5kp05)+k_0*s));
    dE = PHI.*(1+cos(gamma))./s;
    E = sum(sum(dE));
    Intensity(i) = abs(E)^2;
    
end

plot(yi,Intensity);
% 
% %%
% n_wvPts = 12000;
% r_mask = 500;
% l   = linspace(-r_mask,r_mask,n_wvPts);
% lx  = repmat(l,n_wvPts,1);
% ly  = rot90(lx,1);
% lr = sqrt(lx.^2+ly.^2);
% lphi = atan2(ly,lx);
% lphi((n_wvPts/2+1):n_wvPts,:) = lphi((n_wvPts/2+1):n_wvPts,:)+2*pi;
% 
% lPSI_m5kp10 = kp_10*lr+5*lphi;
% lmask_m5kp10 = false(n_wvPts);
% lmask_m5kp10(rem(lPSI_m5kp10,2*pi) > pi) = true; 
% 
% msamp = double(mask_m5kp10);
% new_mask = zeros(n_wvPts);
% t = n_wvPts/n_pts;
% for i = 1:n_pts; sx = msamp(i,:); ssx = interp1(x,sx,l,'nearest'); sb=repmat(ssx,t,1); new_mask((t*(i-1)+1):(t*i),:) = sb; end;
% %imask_m5kp10 = interp2(mask_m5kp10,5);
% 
% %%
% 
% yi = linspace(0,300,301);
% Intensity10000 = zeros(size(yi));
% z = 27000;
% 
% for i = 1:length(yi);
%     display(i);
%     sqr1 = lx.^2+(ly-yi(i)).^2;
%     gamma = atan(sqrt(sqr1)/z);
%     s = sqrt(z^2+sqr1);
%     PHI=exp(1i*(lPSI_m5kp10+k_0*s));
%     %PHI=exp(1i*(pi*double(lmask_m5kp10)+k_0*s));
%     %PHI=exp(1i*(pi*new_mask+k_0*s));
%     dE = PHI.*(1+cos(gamma))./s;
%     E = sum(sum(dE));
%     Intensity10000(i) = abs(E)^2;
%     
% end
% 
% plot(yi,Intensity10000);