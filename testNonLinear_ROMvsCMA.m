close all
clear
setPlot

%% IEC Site definition

Vref      = 60;   % Reference wind speed 
zref      = 87.6; % Hub height above ground
Uref      = 12;   % Wind speed at zref
Iuref     = 0.15; % Turbulence intensity
zmin      = 5;    % Minimum height for wind shear. From 0 to zmin, U is constant and equal to U(z=zmin).

site = siteIEC3(Uref,zref,...
                'category',Iuref,...
                'windTurbineClass', Vref,...
                'zmin',zmin);
%% Definition of time vector 
% Time length of the time series is fixed to be 600 seconds and the
% sampling period is set, according to Shanon theorem, to 1/(2*fmax) being
% fmax the maximum frequency of interest of the PSD. For this particular
% case the maximum frequency is 10 Hz.
tsim    = 60;
fmax    = 10;
dt      = 1/(2*fmax);
t       = 0:dt:tsim-dt;

%% Vertical mesh
nr      = 11; % Number of rows
R       = 63; % Domain semilenght

mesh_v = struct(...
'R',zref,...    
'domainDimension',[0,zref],...
'zmin',0,...
'N_i',[1,nr]);

%% Generation of the wind time series
wstowss_v  = getSPMV1Dts(t,site,mesh_v,...
                         'numberOfRealizations',1,...    
                         'verticalMeshMode','upward');
                     
utilde = wstowss_v.utilde_i;
u      = wstowss_v.u_i;
U      = wstowss_v.U_1;

utilde = permute(utilde,[5,1,2,3,4]);
u      = permute(u,[5,1,2,3,4]);
%% Get the Euler-Bernoulli beam
tower     = nrel5MWTower;
towerBeam = nrel5MWebbTower(tower);
zeta_s = [0.01];

%% Get ROM & modal model
% Define the number of modes and rotational speed
Nmf              = 1;
Omega0           = 0;

% Obtain ROM & modal model
modal  = ebb2modal(towerBeam,Nmf,Omega0,zeta_s);
rom    = ebb2rom(towerBeam,Nmf,Omega0,zeta_s,site);

%% Compute response
ic = zeros(nr,2);
mp = struct('U_1',U);
nA = 5;
analysisVec = linspace(0.8,1.2,nA);

% omega2_ij - Modal model
omegaVec =  analysisVec * modal.eigenValue;
for i = 1:nA 
    omega2Saved = modal.omega2_ij;
    modal.omega2_ij = omegaVec(i)^2;
    
    ebbdr_ad2 = getEbbDynamicResponse(modal,site,t,utilde,ic,mp,...
                                      'modal2modal',@modal2modal_dof1ad2,...
                                      'odemodel',@odelin2forced_dof1ad2);
    modal.omega2_ij = omega2Saved;
    
    wStatsModal  = stowisesNv1d2soStats(ebbdr_ad2.t,ebbdr_ad2.w);
    mbStatsModal = stowisesNv1d2soStats(ebbdr_ad2.t,ebbdr_ad2.mb);

    wMeanModal(i)  = wStatsModal.mean(nr);
    mbMeanModal(i) = mbStatsModal.mean(1);
    wVarModal(i)   = wStatsModal.var(nr);
    mbVarModal(i)  = mbStatsModal.var(1);
    end

% omega2_ij - ROM model non linear
ic = [0 0];

for i = 1:nA 
    omegaSaved   = rom.wT;
    rom.wT = omegaVec(i);
    
    romdr_ad2     = getRomDynamicResponse(rom,site,t,utilde,ic,mp);
    
    rom.wT = omegaSaved;
    
    wStatsRom  = stowisesNv1d2soStats(romdr_ad2{1}.t,romdr_ad2{1}.w);
    mbStatsRom = stowisesNv1d2soStats(romdr_ad2{1}.t,romdr_ad2{1}.mb);

    wMeanRom(i)  = wStatsRom.mean;
    mbMeanRom(i) = mbStatsRom.mean;
    wVarRom(i)   = wStatsRom.var;
    mbVarRom(i)  = mbStatsRom.var;
end

% omega2_ij - ROM model linear
L0_ad = 340/3.87;
% ndRom = rom2ndRom(rom,site);

for i = 1:nA 
    omegaSaved   = rom.wT;
    rom.wT = omegaVec(i);
    ndRom = rom2ndRom(rom,site);
    
    Db_ad = modal.c(1)/modal.c(end);
    alpha = site.alpha;
    Iwt = 0.5 * (1 + alpha + Db_ad/2)/(alpha^2 + 5*alpha +3);
    
    dispMeanAd_Analytic(i) = 0.5* ndRom.gamma * ndRom.lambda^2 * Iwt /(2*pi*ndRom.kapaT)^2;
    bmMeanAd_Analytic(i)   = 0.5 * ndRom.lambda^2 * Iwt;
    
    [dispMeanAd(i),bmMeanAd(i),dispVarAd(i),bmVarAd(i)] = ndRomLinearSolution(ndRom,site,L0_ad);
    
    rom.wT = omegaSaved;
end
kapaTVec = omegaVec .*3.87./(12.*2.*pi);
% Desplazamientos mayores en el modal
% Momentos mayores en el modal
figure(1)
subplot(2,2,1)
plot(kapaTVec,wMeanModal/ndRom.D0,'b-o');hold on
plot(kapaTVec,wMeanRom/ndRom.D0,'r-o');hold on
plot(kapaTVec,dispMeanAd,'k-o');hold on
plot(kapaTVec,dispMeanAd_Analytic,'y-o');
xlabel('$\kappa_{T}$[-]');
ylabel('$\overline{\epsilon_{y}} \Lambda$[-]');
legend('Modal','Rom - Non linear','Rom - Linear','Rom - Linear Analytic')
subplot(2,2,2)
plot(kapaTVec,mbMeanModal/(site.density*ndRom.D0^3*site.Uref^2),'b-o');hold on
plot(kapaTVec,mbMeanRom/(site.density*ndRom.D0^3*site.Uref^2),'r-o');hold on
plot(kapaTVec,bmMeanAd,'k-o');
plot(kapaTVec,bmMeanAd_Analytic,'y-o');
xlabel('$\kappa_{T}$[-]');
ylabel('$M_{y_{G}}^{r,r}/ M_{C}U_{0}^{2}$[-]');
legend('Modal','Rom - Non linear','Rom - Linear','Rom - Linear Analytic')
subplot(2,2,3) 
plot(kapaTVec,wVarModal/ndRom.D0^2,'b-o');hold on
plot(kapaTVec,wVarRom/ndRom.D0^2,'r-o');hold on
plot(kapaTVec,dispVarAd,'k-o');
xlabel('$\kappa_{T}$[-]');
ylabel('$\sigma_{\epsilon_{y}H}^{2} /D_{0}^{2}$[-]');
legend('Modal','Rom - Non linear','Rom - Linear')
subplot(2,2,4)
plot(kapaTVec,mbVarModal/(site.density*ndRom.D0^3*site.Uref^2)^2,'b-o');hold on
plot(kapaTVec,mbVarRom/(site.density*ndRom.D0^3*site.Uref^2)^2,'r-o');hold on
plot(kapaTVec,bmVarAd,'k-o');
xlabel('$\kappa_{T}$[-]');
ylabel('$\sigma^{2}_{m_{y_{G}}^{r,r}}/ M_{C}^{2}U_{0}^{4}$[-]');
legend('Modal','Rom - Non linear','Rom - Linear')

figure(11) % Displacement
plot(ebbdr_ad2.t,ebbdr_ad2.w(:,nr),'b-'); hold on;
plot(romdr_ad2{1}.t,romdr_ad2{1}.w,'r-'); %
xlabel('$t$ [s]'); 
ylabel('$w$ [m]');
legend('Modal','Rom')

figure(12) % Velocity
plot(ebbdr_ad2.t,ebbdr_ad2.dwdt(:,nr),'b-'); hold on;
plot(romdr_ad2{1}.t,romdr_ad2{1}.dwdt,'r-'); %
xlabel('$t$ [s]'); 
ylabel('$\dot{w} [m/s]$');
legend('Modal','Rom')

figure(13) % Acceleration
plot(ebbdr_ad2.t,ebbdr_ad2.d2wdt2(:,nr),'b-'); hold on;
plot(romdr_ad2{1}.t,romdr_ad2{1}.d2wdt2,'r-'); 
xlabel('$t$ [s]'); 
ylabel('$\ddot{w}[m/s^{2}]$');
legend('Modal','Rom')

figure(14) % Wind force
plot(ebbdr_ad2.t,ebbdr_ad2.fj,'b-'); hold on;
plot(romdr_ad2{1}.t,romdr_ad2{1}.fj,'r-'); 
xlabel('$t$ [s]'); 
ylabel('$f_{z_{G}}$[N]');
legend('Modal','Rom')

figure(15) % Bending moment
plot(ebbdr_ad2.t,ebbdr_ad2.mb(:,1),'b-'); hold on;
plot(romdr_ad2{1}.t,romdr_ad2{1}.mb,'r-'); 
xlabel('$t$ [s]'); 
ylabel('$M_{y_{G}}^{r,r}$[Nm]');
legend('Modal','Rom')

% Plot PSDs
nt          = length(t);

[Sq_m,fq_m]  = pwelch(ebbdr_ad2.w(:,nr),rectwin(nt),0,nt,1/dt);
[Sv_m,fv_m]  = pwelch(ebbdr_ad2.dwdt(:,nr),rectwin(nt),0,nt,1/dt);
[Sa_m,fa_m]  = pwelch(ebbdr_ad2.d2wdt2(:,nr),rectwin(nt),0,nt,1/dt);
[Sfj_m,ffj_m] = pwelch(ebbdr_ad2.fj,rectwin(nt),0,nt,1/dt);
[Smb_m,fmb_m] = pwelch(ebbdr_ad2.mb(:,1),rectwin(nt),0,nt,1/dt);

[Sq_r,fq_r]  = pwelch(romdr_ad2{1}.w,rectwin(nt),0,nt,1/dt);
[Sv_r,fv_r]  = pwelch(romdr_ad2{1}.dwdt,rectwin(nt),0,nt,1/dt);
[Sa_r,fa_r]  = pwelch(romdr_ad2{1}.d2wdt2,rectwin(nt),0,nt,1/dt);
[Sfj_r,ffj_r] = pwelch(romdr_ad2{1}.fj,rectwin(nt),0,nt,1/dt);
[Smb_r,fmb_r] = pwelch(romdr_ad2{1}.mb,rectwin(nt),0,nt,1/dt);

figure(21) % Displacement PSD
loglog(fq_m,fq_m.*Sq_m,'b-'); hold on; 
loglog(fq_r,fq_r.*Sq_r,'r-');
xlabel('$f$ [Hz]'); 
ylabel('$fS_{w}[m^{2}]$');
grid on
legend('Modal','Rom')

figure(22) % Velocity PSD
loglog(fv_m,fv_m.*Sv_m,'b-'); hold on; 
loglog(fv_r,fv_r.*Sv_r,'r-');
xlabel('$f$ [Hz]'); 
ylabel('$fS_{\dot{w}} [m^{2}/s^{2}]$');
grid on
legend('Modal','Rom')

figure(23) % Acceleration PSD
loglog(fa_m,fa_m.*Sa_m,'b-'); hold on; 
loglog(fa_r,fa_r.*Sa_r,'r-');
xlabel('$f$ [Hz]'); 
ylabel('$fS_{\ddot{w}} [m^{2}/s^{4}]$');
grid on
legend('Modal','Rom')

figure(24) % Wind force PSD
loglog(ffj_m,ffj_m.*Sfj_m,'b-'); hold on; 
loglog(ffj_r,ffj_r.*Sfj_r,'r-');
xlabel('$f$ [Hz]'); 
ylabel('$fS_{f_{z_{G}}} [N^{2}]$');
grid on
legend('Modal','Rom')

figure(25) % Wind force PSD
loglog(fmb_m,fmb_m.*Smb_m,'b-'); hold on; 
loglog(fmb_r,fmb_r.*Smb_r,'r-');
xlabel('$f$ [Hz]'); 
ylabel('$fS_{M_{y_{G}}^{r,r}} [N^{2}m^{2}]$');
grid on
legend('Modal','Rom')

