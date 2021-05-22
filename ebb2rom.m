function rom = ebb2rom(ebb,Nmf,Omega0,zeta_s,site)
%EBB2MODAL  Computes ROM parameters of a blade of a wind turbine
%
%   ROM = EBB2ROM(B,N,W,Z) computes ROM parameters structure ROM of
%   a 1D Euler-Bernoulli beam, B, with angular velocity W corresponding
%   to the first N modes with a Nx1 column vector of structural modal 
%   damping ratios Z. The output structure of ROM parameters, ROM, has 
%   the following fields:
%           mT: double
%           xGT: double
%           IT: double
%           wT: double
%    omega2_ij: double
%           kT: double
%           cT: double
%           cs: double
%            r: [nx1 double]
%            D: [nx1 double]
%           kY: double
%         zeta: double
%       thetaG: [nx1 double]
%      airfoil: [nx1 double]

m  = ebb2modal(ebb,Nmf,Omega0,zeta_s);

% Mass
mT  = trapz(ebb.r,ebb.m);

% Lengths
xGT = trapz(ebb.r,ebb.m .* ebb.r)/mT;
r  = ebb.r;
D =  ebb.c;

% Inertia
IT  = trapz(ebb.r,ebb.m.*ebb.r.^2);

% Stiffness
wT        = m.eigenValue(1,1);
omega2_ij = wT^2;
kT        = wT^2 * IT;
kY        = site.gravity * xGT * mT + kT;

%Damping
zeta = m.zeta_s(1,1);
cs   = 2 * zeta * wT;
cT   = 2 * zeta * IT * wT;

% Airfoil
thetaG  = ebb.thetaG;
airfoil = ebb.airfoil;

rom     = struct('mT',mT,...
                   'xGT',xGT,...
                   'IT',IT,...
                   'wT',wT,...  
                   ...'omega2_ij',omega2_ij,... 
                   ...'kT',kT,...  
                   ...'cT',cT,...
                   ...'cs',cs,...
                   'r',r,...
                   'D',D,...
                   ...'kY',kY,...
                   'zeta',zeta,...
                   'thetaG',thetaG);   
               
rom.airfoil = ebb.airfoil;
end