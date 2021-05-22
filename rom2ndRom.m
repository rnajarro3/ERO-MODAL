function romAd = rom2ndRom(rom,site)
% ROM2NDROM  Computes non dimensional ROM parameters of a blade of a wind turbine
%
%   ROMAD = ROM2NDROM(ROM,U,SITE) computes non dimensional ROM parameters structure of
%   a ROM, with mean wind velocity U and site SITE. 
%   The output structure of ROM parameters, ROM, has the following fields:
%           t_c: double
%           M_C: double
%         mT_ad: double
%        lambda: double
%            D0: double
%        xGT_ad: double
%          x_ad: [nx1 double]
%          D_ad: [nx1 double]
%         Lu_ad: double
%         gamma: double
%         wT_ad: double
%         kapaT: double
%         kY_ad: double
%         cT_ad: double
%          zeta: double
%            CD: double

H   = site.zref;
rho = site.density;
U0  = site.Uref;
D   = rom.D;
r   = rom.r;

D0 = interp1(r,D,H);

[nr,nc] = size(D);
alpha             = pi/2 - rom.thetaG;
cd_i              = zeros(nr,1);
for i = 1:nr
    fcd                = rom.airfoil{i}.cd;
    cd_i(i)            = fcd(alpha(i));
end

% Time
t_c    = D0 /U0;

%Mass 
M_C    = rho * D0^3;
mT_ad  = rom.mT / M_C;

% Lenghts
lambda = H/D0;
xGT_ad = rom.xGT / H;
D_ad   = D/D0;
x_ad   = r / H;

% Lock
gamma  = rho * H * D0^4/ rom.IT;

% Stiffness
kapaT  = rom.wT*D0/ U0/(2*pi);

%Damping
zeta   = rom.zeta;

romAd     = struct(...'t_c',t_c,...
                   'D0',D0,...
                   'mT_ad',mT_ad,...
                   'lambda',lambda,...  
                   'xGT_ad',xGT_ad,...
                   'x_ad',x_ad,...
                   'D_ad',D_ad,...
                   ...'Lu_ad',Lu_ad,...
                   'gamma',gamma,...
                   ...'wT_ad',wT_ad,...
                   'kapaT',kapaT,...
                   ...'kY_ad',kY_ad,...
                   ...'cT_ad',cT_ad,...
                   'zeta',zeta,...
                   'CD',cd_i);   
end