close all
clear
setPlot
line    = {'k-o','r-o','b-o','m-o','g-o',...
           'k--o','r--o','b--o','m--o','g--o'};
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
tsim    = 600;
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

u      = wstowss_v.u_i;
U      = wstowss_v.U_1;

u = permute(u,[5,1,2,3,4]);
%% Get the Euler-Bernoulli beam
tower     = nrel5MWTower;
towerBeam = nrel5MWebbTower(tower);
zeta_s = [0.01];

%% Get ROM 
% Define the number of modes and rotational speed
Nmf              = 1;
Omega0           = 0;
% Obtain ROM
rom = ebb2rom(towerBeam,Nmf,Omega0,zeta_s,site);
% Obtain ndROM
ndRom = rom2ndRom(rom,site);

%% 2D PLOTS
ic = [0 0];
Na = 1;
L0_ad = 340/3.87;
nA = 5;
analysisVec  = linspace(0.8,1.2,nA);
analysisVec2 = linspace(0.95,1.05,nA);

% Alpha
alphaVec =  analysisVec * site.alpha;
for i = 1:nA 
    alphaSaved = site.alpha;
    site.alpha = alphaVec(i);
    
    [dispMeanAd(i),bmMeanAd(i),dispVarAd(i),bmVarAd(i),Ieps(i),IdispVar(i),IbmVar(i)] = ndRomLinearSolution(ndRom,site,L0_ad);
    
    site.alpha = alphaSaved;
end
xAxis = '$\alpha$ [-]';
plotMean(11,alphaVec,dispMeanAd,bmMeanAd,IdispVar,xAxis)
% plotVariance(2,alphaVec,dispVarAd,bmVarAd,IdispVar,IbmVar,xAxis)
% 
% % lambda
% lambdaVec = analysisVec * 22;
% for i = 1:nA 
%     ndRom = rom2ndRom(rom,site);
%     lambdaSaved = ndRom.lambda;
%     ndRom.lambda = lambdaVec(i);
%     [dispMeanAd(i),bmMeanAd(i),dispVarAd(i),bmVarAd(i),Ieps(i),IdispVar(i),IbmVar(i)] = ndRomLinearSolution(ndRom,site,L0_ad);
%     ndRom.lambda = lambdaSaved;
% end
% xAxis = '$\Lambda$ [-]';
% plotVariance(3,lambdaVec,dispVarAd,bmVarAd,IdispVar,IbmVar,xAxis)
% 
% % Lu/D0
% L0Vec = analysisVec * L0_ad;
% for i = 1:nA 
%     ndRom = rom2ndRom(rom,site);
%     L0Saved = L0_ad;
%     L0_ad = L0Vec(i);
%     [dispMeanAd(i),bmMeanAd(i),dispVarAd(i),bmVarAd(i),Ieps(i),IdispVar(i),IbmVar(i)] = ndRomLinearSolution(ndRom,site,L0_ad);
%     L0_ad = L0Saved;
% end
% xAxis = '$L/D_{0}$ [-]';
% plotVariance(4,L0Vec,dispVarAd,bmVarAd,IdispVar,IbmVar,xAxis)
% 
% % kappaT
% wTVec = analysisVec * rom.wT;
% for i = 1:nA 
%     wTSaved = rom.wT;
%     rom.wT = wTVec(i);
%     ndRom = rom2ndRom(rom,site);
%     [dispMeanAd(i),bmMeanAd(i),dispVarAd(i),bmVarAd(i),Ieps(i),IdispVar(i),IbmVar(i)] = ndRomLinearSolution(ndRom,site,L0_ad);
%     rom.wT = wTSaved;
%     kapaTVec(i) = ndRom.kapaT;
% end
% xAxis = '$\kappa_{T}$ [-]';
% plotVariance(5,kapaTVec,dispVarAd,bmVarAd,IdispVar,IbmVar,xAxis)
% 
% % zetaS
% zetaVec = analysisVec * rom.zeta;
% for i = 1:nA 
%     zetaSaved = rom.zeta;
%     rom.zeta = zetaVec(i);
%     ndRom = rom2ndRom(rom,site);
%     [dispMeanAd(i),bmMeanAd(i),dispVarAd(i),bmVarAd(i),Ieps(i),IdispVar(i),IbmVar(i)] = ndRomLinearSolution(ndRom,site,L0_ad);
%     rom.zeta = zetaSaved;
% end
% xAxis = '$\zeta_{s}$ [-]';
% plotVariance(6,zetaVec,dispVarAd,bmVarAd,IdispVar,IbmVar,xAxis)

%% 3D PLOTS
%1. lambda - alpha
alphaVec  = analysisVec * site.alpha;
lambdaVec = analysisVec2 * ndRom.lambda;
for i = 1:nA 
    for j = 1:nA
    alphaSaved = site.alpha;
    site.alpha = alphaVec(j);
    
    lambdaSaved = ndRom.lambda;
    ndRom.lambda = lambdaVec(i);
    
    [dispMeanAd(i,j),bmMeanAd(i,j),dispVarAd(i,j),bmVarAd(i,j),Ieps(i,j),IdispVar(i,j),IbmVar(i,j)] = ndRomLinearSolution(ndRom,site,L0_ad);
    
    site.alpha = alphaSaved;
    ndRom.lambda = lambdaSaved;
    end
end
xAxis = '$\Lambda$ [-]';
leg   = '$\alpha$';
plotVarianceParam(1,lambdaVec,alphaVec,dispVarAd,bmVarAd,IdispVar,IbmVar,xAxis,leg,line)

%2. lambda - L0/DO
lambdaVec = analysisVec2 * ndRom.lambda;
L0Vec    = analysisVec * L0_ad;
for i = 1:nA 
    for j = 1:nA
    lambdaSaved = ndRom.lambda;
    ndRom.lambda = lambdaVec(i);
    
    L0Saved = L0_ad;
    L0_ad = L0Vec(j);
    
    [dispMeanAd(i,j),bmMeanAd(i,j),dispVarAd(i,j),bmVarAd(i,j),Ieps(i,j),IdispVar(i,j),IbmVar(i,j)] = ndRomLinearSolution(ndRom,site,L0_ad);
    
    ndRom.lambda = lambdaSaved;
    L0_ad = L0Saved;
    end
end
xAxis = '$\Lambda$ [-]';
leg   = '$L_{0} / D_{0}$ [-]';
plotVarianceParam(2,lambdaVec,L0Vec,dispVarAd,bmVarAd,IdispVar,IbmVar,xAxis,leg,line)

%3. lambda - kapaT
lambdaVec = analysisVec2 * ndRom.lambda;
kapaTVec  = analysisVec * ndRom.kapaT;
for i = 1:nA 
    for j = 1:nA
    lambdaSaved = ndRom.lambda;
    ndRom.lambda = lambdaVec(i);
    
    kapaTSaved = ndRom.kapaT;
    ndRom.kapaT = kapaTVec(j);
    
    [dispMeanAd(i,j),bmMeanAd(i,j),dispVarAd(i,j),bmVarAd(i,j),Ieps(i,j),IdispVar(i,j),IbmVar(i,j)] = ndRomLinearSolution(ndRom,site,L0_ad);
    
    ndRom.lambda = lambdaSaved;
    ndRom.kapaT = kapaTSaved;
    end
end
xAxis = '$\Lambda$ [-]';
leg   = '$\kappa_{T}$';
plotVarianceParam(3,lambdaVec,kapaTVec,dispVarAd,bmVarAd,IdispVar,IbmVar,xAxis,leg,line)

%4. lambda - zetaS
lambdaVec = analysisVec2 * ndRom.lambda;
zetaVec    = analysisVec * ndRom.zeta;
for i = 1:nA 
    for j = 1:nA
    lambdaSaved = ndRom.lambda;
    ndRom.lambda = lambdaVec(i);
    
    zetaSaved = ndRom.zeta;
    ndRom.zeta = zetaVec(j);
    
    [dispMeanAd(i,j),bmMeanAd(i,j),dispVarAd(i,j),bmVarAd(i,j),Ieps(i,j),IdispVar(i,j),IbmVar(i,j)] = ndRomLinearSolution(ndRom,site,L0_ad);
    
    ndRom.lambda = lambdaSaved;
    ndRom.zeta = zetaSaved;
    end
end
xAxis = '$\Lambda$ [-]';
leg   = '$\zeta_{S}$';
plotVarianceParam(4,lambdaVec,zetaVec,dispVarAd,bmVarAd,IdispVar,IbmVar,xAxis,leg,line)

%5. Alpha - L0/D0
alphaVec = analysisVec * site.alpha;
L0Vec    = analysisVec * L0_ad;
for i = 1:nA 
    for j = 1:nA
    alphaSaved = site.alpha;
    site.alpha = alphaVec(i);
    
    L0Saved = L0_ad;
    L0_ad = L0Vec(j);
    
    [dispMeanAd(i,j),bmMeanAd(i,j),dispVarAd(i,j),bmVarAd(i,j),Ieps(i,j),IdispVar(i,j),IbmVar(i,j)] = ndRomLinearSolution(ndRom,site,L0_ad);
    
    site.alpha = alphaSaved;
    L0_ad = L0Saved;
    end
end
xAxis = '$\alpha$ [-]';
leg   = '$L_{0} / D_{0}$ [-]';
plotVarianceParam(5,alphaVec,L0Vec,dispVarAd,bmVarAd,IdispVar,IbmVar,xAxis,leg,line)

%6. Alpha - kapaT
alphaVec = analysisVec * site.alpha;
kapaTVec   = analysisVec * ndRom.kapaT;
for i = 1:nA 
    for j = 1:nA
    alphaSaved = site.alpha;
    site.alpha = alphaVec(i);
    
    kapaTSaved = ndRom.kapaT;
    ndRom.kapaT = kapaTVec(j);
    
    [dispMeanAd(i,j),bmMeanAd(i,j),dispVarAd(i,j),bmVarAd(i,j),Ieps(i,j),IdispVar(i,j),IbmVar(i,j)] = ndRomLinearSolution(ndRom,site,L0_ad);
    
    site.alpha = alphaSaved;
    ndRom.kapaT = kapaTSaved;
    end
end
xAxis = '$\alpha$ [-]';
leg   = '$\kappa_{T}$';
plotVarianceParam(6,alphaVec,kapaTVec,dispVarAd,bmVarAd,IdispVar,IbmVar,xAxis,leg,line)

%7. Alpha - zetaS
alphaVec = analysisVec * site.alpha;
zetaVec    = analysisVec * ndRom.zeta;
for i = 1:nA 
    for j = 1:nA
    alphaSaved = site.alpha;
    site.alpha = alphaVec(i);
    
    zetaSaved = ndRom.zeta;
    ndRom.zeta = zetaVec(j);
    
    [dispMeanAd(i,j),bmMeanAd(i,j),dispVarAd(i,j),bmVarAd(i,j),Ieps(i,j),IdispVar(i,j),IbmVar(i,j)] = ndRomLinearSolution(ndRom,site,L0_ad);
    
    site.alpha = alphaSaved;
    ndRom.zeta = zetaSaved;
    end
end
xAxis = '$\alpha$ [-]';
leg   = '$\zeta_{S}$';
plotVarianceParam(7,alphaVec,zetaVec,dispVarAd,bmVarAd,IdispVar,IbmVar,xAxis,leg,line)

%8. Lu/D0 - kapaT
L0Vec    = analysisVec * L0_ad;
kapaTVec = analysisVec * ndRom.kapaT;

for i = 1:nA 
    for j = 1:nA
    L0Saved = L0_ad;
    L0_ad = L0Vec(i);
    
    kapaTSaved = ndRom.kapaT;
    ndRom.kapaT = kapaTVec(j);
    
    [dispMeanAd(i,j),bmMeanAd(i,j),dispVarAd(i,j),bmVarAd(i,j),Ieps(i,j),IdispVar(i,j),IbmVar(i,j)] = ndRomLinearSolution(ndRom,site,L0_ad);
    
    L0_ad = L0Saved;
    ndRom.kapaT = kapaTSaved;
    end
end
xAxis = '$L_{0} / D_{0}$ [-]';
leg   = '$\kappa_{T}$ [-]';
plotVarianceParam(8,L0Vec,kapaTVec,dispVarAd,bmVarAd,IdispVar,IbmVar,xAxis,leg,line)

%9. Lu/D0 - zetaS
L0Vec   = analysisVec * L0_ad;
zetaVec = analysisVec * ndRom.zeta;

for i = 1:nA 
    for j = 1:nA
    L0Saved = L0_ad;
    L0_ad = L0Vec(i);
    
    zetaSaved = ndRom.zeta;
    ndRom.zeta = zetaVec(j);
    
    [dispMeanAd(i,j),bmMeanAd(i,j),dispVarAd(i,j),bmVarAd(i,j),Ieps(i,j),IdispVar(i,j),IbmVar(i,j)] = ndRomLinearSolution(ndRom,site,L0_ad);
    
    L0_ad = L0Saved;
    ndRom.zeta = zetaSaved;
    end
end
xAxis = '$L_{0} / D_{0}$ [-]';
leg   = '$\zeta_{S}$ [-]';
plotVarianceParam(9,L0Vec,zetaVec,dispVarAd,bmVarAd,IdispVar,IbmVar,xAxis,leg,line)

%10. kapaT - zetaS
kapaTVec = analysisVec * ndRom.kapaT;
zetaVec  = analysisVec * ndRom.zeta;

for i = 1:nA 
    for j = 1:nA
    kapaTSaved = ndRom.kapaT;
    ndRom.kapaT = kapaTVec(i);
    
    zetaSaved  = ndRom.zeta;
    ndRom.zeta = zetaVec(j);
    
    [dispMeanAd(i,j),bmMeanAd(i,j),dispVarAd(i,j),bmVarAd(i,j),Ieps(i,j),IdispVar(i,j),IbmVar(i,j)] = ndRomLinearSolution(ndRom,site,L0_ad);
    
    ndRom.kapaT = kapaTSaved;
    ndRom.zeta  = zetaSaved;
    end
end
xAxis = '$\kappa_{T}$ [-]';
leg   = '$\zeta_{S}$ [-]';
plotVarianceParam(10,kapaTVec,zetaVec,dispVarAd,bmVarAd,IdispVar,IbmVar,xAxis,leg,line)