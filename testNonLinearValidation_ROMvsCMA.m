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

modal.m_ij      = rom.IT/rom.r(end)^2;
modal.Inv_m_ij  = rom.r(end)/rom.IT;
modal.eigenMode = rom.r/rom.r(end);
modal.omega2_ij = rom.omega2_ij;
modal.cs_ij     = rom.cs;

%% Compute response
%Modal model
ic = zeros(nr,2);
mp = struct('U_1',U);
ebbdr_ad2 = getEbbDynamicResponse(modal,site,t,utilde,ic,mp,...
                                  'modal2modal',@modal2modal_dof1ad2,...
                                  'odemodel',@odelin2forced_dof1ad2);

% ROM model
ic = [0 0];
romdr_ad2 = getRomDynamicResponse(rom,site,t,utilde,ic,mp);
    
%% Plots
figure(11) % Displacement
plot(ebbdr_ad2.t,ebbdr_ad2.w(:,nr)*rom.r(end),'b-'); hold on;
plot(romdr_ad2{1}.t,romdr_ad2{1}.w,'r-');
xlabel('$t$ [s]'); 
ylabel('$w$ [m]');
legend('Modal','Rom')

figure(12) % Velocity
plot(ebbdr_ad2.t,ebbdr_ad2.dwdt(:,nr)*rom.r(end),'b-'); hold on;
plot(romdr_ad2{1}.t,romdr_ad2{1}.dwdt,'r-');
xlabel('$t$ [s]'); 
ylabel('$\textrm{d} w/ \textrm{d} t(t)$[m/s]');
legend('Modal','Rom')

figure(13) % Acceleration
plot(ebbdr_ad2.t,ebbdr_ad2.d2wdt2(:,nr)*rom.r(end),'b-'); hold on;
plot(romdr_ad2{1}.t,romdr_ad2{1}.d2wdt2,'r-'); 
xlabel('$t$ [s]'); 
ylabel('$\textrm{d}^{2} w/ \textrm{d} t^{2}(t)$[m/s2]');
legend('Modal','Rom')

figure(14) % Wind force
plot(ebbdr_ad2.t,ebbdr_ad2.fj,'b-'); hold on;
plot(romdr_ad2{1}.t,romdr_ad2{1}.fj,'r-'); 
xlabel('$t$ [s]'); 
ylabel('$ft(t)$[N]');
legend('Modal','Rom')

figure(15) % Bending moment
plot(ebbdr_ad2.t,ebbdr_ad2.mb(:,1)*rom.r(end),'b-'); hold on;
plot(romdr_ad2{1}.t,romdr_ad2{1}.mb,'r-'); 
xlabel('$t$ [s]'); 
ylabel('$mb(t)$[Nm]');
legend('Modal','Rom')

%% Plot PSDs
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
loglog(fq_m,fq_m.*Sq_m*rom.r(end)^2,'b-'); hold on; 
loglog(fq_r,fq_r.*Sq_r,'r-');
legend('Modal','Rom')

figure(22) % Velocity PSD
loglog(fv_m,fv_m.*Sv_m*rom.r(end)^2,'b-'); hold on; 
loglog(fv_r,fv_r.*Sv_r,'r-');
legend('Modal','Rom')

figure(23) % Acceleration PSD
loglog(fa_m,fa_m.*Sa_m*rom.r(end)^2,'b-'); hold on; 
loglog(fa_r,fa_r.*Sa_r,'r-');
legend('Modal','Rom')

figure(24) % Wind force PSD
loglog(ffj_m,ffj_m.*Sfj_m,'b-'); hold on; 
loglog(ffj_r,ffj_r.*Sfj_r,'r-');
legend('Modal','Rom')

figure(25) % Bending moment PSD
loglog(fmb_m,fmb_m.*Smb_m*rom.r(end)^2,'b-'); hold on; 
loglog(fmb_r,fmb_r.*Smb_r,'r-');
legend('Modal','Rom')

