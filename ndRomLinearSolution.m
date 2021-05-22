function [dispMeanAd,bmMeanAd,dispVarAd,bmVarAd,Ieps,IdispVar,IbmVar]  = ndRomLinearSolution(ndRom,site,L0_ad)

U0   = site.Uref;
H    = site.zref;
Iu   = site.Iuref;

D0     = ndRom.D0;
zeta_s = ndRom.zeta;
lambda = ndRom.lambda;
gamma  = ndRom.gamma;
kapaT  = ndRom.kapaT;
D_ad   = ndRom.D_ad;
x_ad   = ndRom.x_ad;
CD     = ndRom.CD;
xGT_ad = ndRom.xGT_ad;
mT_ad  = ndRom.mT_ad;

% Compute non dimensional parameters and variables
t_c   = D0/U0;
wT_ad = 2*pi*kapaT;
kY_ad = wT_ad^2 *lambda/gamma  + (site.gravity*H/U0^2) * xGT_ad * mT_ad;
cY_ad = 2 * zeta_s * wT_ad *lambda/gamma;
U_ad  = powerWindShear(x_ad*H,U0,H,site.alpha,5)/U0;

% Mean top displacement and mean bending moment
Ieps       = trapz(x_ad,CD .* D_ad .* U_ad .^ 2 .* x_ad );
dispMeanAd = gamma / wT_ad^2 * lambda^2 / 2 * Ieps;
bmMeanAd   =                   lambda^2 / 2 * Ieps;

% Variance of top displacement and bending moment
fmin = log10(1e-4);
fmax = log10(10);
nf   = 501;
f_i  = logspace(fmin,fmax,nf)';
kapa = t_c * f_i;

zeta_a = gamma / wT_ad * lambda^2 / 2 * trapz(x_ad,CD .* D_ad .* U_ad .^ 2 .* x_ad .^2 );
zeta_T = zeta_a + zeta_s;
HT_ad  = gamma * lambda./abs(-kapa.^2 + 2*sqrt(-1)*zeta_T*kapa.*kapaT + kapaT^2)/(2*pi)^2;

Su_ad  = @(r1,r2,kapa)ndVonKarmanCrossSpectra(kapa*L0_ad,lambda/L0_ad*abs(r1-r2));
% Su_ad     = @(kapa)kaimalSpectra(kapa*L0_ad,1,1,1);
% Coh_ad    = @(r1,r2,kapa)davenportCoherence(lambda/L0_ad*abs(r1-r2),kapa*L0_ad, 1, 1);
% S_f1f1_ad = getSf1f1(kapa,x_ad,x_ad,CD.*D_ad.*U_ad,Su_ad,Coh_ad);
S_f1f1_ad = getSf1f1vK(kapa,x_ad,x_ad,CD.*D_ad.*U_ad,Su_ad);

IdispVar  = trapz(kapa,                                    HT_ad.^2.*S_f1f1_ad);
IbmVar    = trapz(kapa,((2*pi.*kapa.*cY_ad).^2 + kY_ad^2).*HT_ad.^2.*S_f1f1_ad);
dispVarAd = lambda^2 * Iu^2 * L0_ad * IdispVar;
bmVarAd   =            Iu^2 * L0_ad * IbmVar;

end