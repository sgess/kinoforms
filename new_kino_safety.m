gamma_deg_300 = 0.312;
gamma_deg_400 = 0.273;
gamma_deg_500 = 0.250;
gamma_deg_600 = 0.235;

gamma_rad_300 = pi*gamma_deg_300/180;
gamma_rad_400 = pi*gamma_deg_400/180;
gamma_rad_500 = pi*gamma_deg_500/180;
gamma_rad_600 = pi*gamma_deg_600/180;

L = 3987;
off_300 = 17;
off_400 = 14.5;
off_500 = 13;
off_600 = 12;

id_300 = L*tan(gamma_rad_300)-off_300
id_400 = L*tan(gamma_rad_400)-off_400
id_500 = L*tan(gamma_rad_500)-off_500
id_600 = L*tan(gamma_rad_600)-off_600

off = 19;
gamma_rad_ax = (1.45-1)*pi*0.75/180;
id_ax = L*tan(gamma_rad_ax)-off