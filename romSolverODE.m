function roma = romSolverODE(tj,rom,y0)
% ROMA = romSolverODE(TSPAN,ROM,IC) solves the dynamic response of a 1D
% continuous system defined by the ROM data structure during the time 
% defined by TSPAN and with initial conditions IC using
% Reduce Order Model Analysis (ROMA). The output structure ROMA is defined
% as follows
%      t: []
%     fj: []
%      q: []
%      v: []
%      a: []
%  d0eps: []
%  d1eps: []
%  d2eps: []


%% ODE function handle and solution 
odefun        = @(t,y) odeRomModel(t,y,rom);
sol           = ode45(odefun,tj,y0);

%% Evaluation of ode solution
[y,dy]      = mydeval(tj,sol);

d0eps       = y(2,:);
d1eps       = y(1,:);
d2eps       = dy(1,:);

% Transform principal coordinates into correct dimensions
d0eps       = d0eps';
d1eps       = d1eps';
d2eps       = d2eps';

% Get time series of displacements, velocities and accelerations
% rom.Psi_ij = rom.Psi_ij;

q         = d0eps*rom.r(end);
v         = d1eps*rom.r(end);
a         = d2eps*rom.r(end);

mb = rom.IT * rom.wT^2 * q/rom.r(end) +  2 * rom.zeta * rom.IT * rom.wT * v/rom.r(end);
roma       = struct('t',tj,...
                    'fj',rom.fj,...
                    'w',q,...
                    'dwdt',v,...
                    'd2wdt2',a,...
                    'eps',d0eps,...
                    'depsdt',d1eps,...
                    'd2epsdt2',d2eps, ...
                    'mb',mb);

end