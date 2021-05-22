function [Sf1f1,Sf,Su] = getSf1f1vK(fspan,zspan,psi_1,c_a1,spectra)
%getSf1f1  Conputes the PSD of first modal forces
%
%   Sf1f1 = GETSF1F1(FSPAN,RSPAN,PSI,CA,SU,COH) computes 
%   the PSD of first modal forces vector at the frequencies defined by the
%   FSPAN vector of frequencies in Hz, given the first eigenmode vector,
%   PSI, the aerodynamic damping vector, CA, which are functions 
%   on the structural domain defined by RSPAN, the wind PSD function handle
%   SU, and the coherence function handle, COH. The SU function
%   handle is a function of one variable, and it should return a vector of
%   size the same size of the input. The COH function handle is a function
%   of three variables, and it should return a 3D-matrix of size the size
%   of the input arguments. The PSD of first modal forces is the following
%   integral, I:
%
%   I = int_0^R int_0^R ca(r)Psi(r)Su(f)Coh(r,r',f)ca(r')Psi(r') dr dr'
%
%
%   Example of usage: Compute the PSD of the first modal forces for a white
%   noise of intensity 1 and a perfectly correlated wind, i.e. Coh=1
%
%   s0         = 1.0;
%   spectra    = @(f) s0.*ones(size(f));
%   coherence  = @(r1,r2,f) onesCoh(abs(r1-r2),f);
%   psi1       = linspace(0,1,51)';
%   rspan      = linspace(0,1,51)';
%   c_a        = ones(51,1);
%   fspan      = linspace(1e-3,1,1001);
%   Sf1f1      = getSf1f1(fspan,rspan,psi1,c_a,spectra,coherence);
%
%   where onesCoh is a function like this
%   function coh  = onesCoh(dr,f)
%   coh           = ones(size(dr));
%
%   As it can be easily checked, for this very simple example the exact 
%   value of the integral is equal to 1/4.
%
%   A more realistic example is the following one in wich a kaimal spectra
%   and a Davenport coherence function are used. First we define wind
%   speed, density, turbulence intensity and then load aerodynamic and
%   structural data of the NREL-5MW
%
%     U1               = 12;
%     rho              = 1.225;
%     Iu               = 0.15;
%     Lu               = 100;
%     cd               = 0.05;
%     nrel5mwBeam      = beamNREL_5MW49tp;
%     nrel5mwBlade     = @bladeNREL_5MW49tp;
%     [xR,x]           = nrel5mwBlade('xvec');
%     [xR,cadvec]      = nrel5mwBlade('cadvec');
%     chord            = cadvec*R;
%     R                = nrel5mwBeam.radius(end);
%     radius           = x.*R;
%     Nm               = 1;
%     Omega0           = 0;
%     [omega_j,psi_1,d1p_j,d2psi_j] = solveEigen(nrel5mwBeam,Nm,Omega0);
%     fmin             = log10(1e-4);
%     fmax             = log10(10);
%     nf               = 101;
%     f_i              = logspace(fmin,fmax,nf);
%     c_a              = rho.*chord.*cd.*U1;
%     Su               = @(f)kaimalSpectra(f,U1,Lu,Iu*U1);
%     Coh              = @(r1,r2,f)davenportCoherence(abs(r1-r2),f,U1,Lu);
%     S_f1f1           = getSf1f1(f_i,radius,psi_1',c_a,Su,Coh);
%     semilogx(f_i,S_f1f1)

% Get number of angular frecuencies FIXME
nw      = numel(fspan);

% z row vector
% eps is added in order to avoid the singularity of the coh
% function when z = 0
z       = zspan + eps;

% Get number of structural points
nz      = length(zspan);


% ensure c_a row vector
c_a1    = c_a1(:)';

% column vector
c_a2    = c_a1';
c_az1   = c_a1(ones(nz,1),:,ones(nw,1));
c_az2   = c_a2(:,ones(nz,1),ones(nw,1));


% Wind CPSD
Su       = site2cpsd(fspan,z,spectra);

% Force CPSD
Sf       = c_az1.*c_az2.*Su;

% Double structural domain integration of the force CPSD 
Sf1f1    = int2sf(zspan,psi_1,Sf);


function Su = site2cpsd(f,z,spectra)
% PSD is one-sided and frequency dependent

% Get number of angular frecuencies
nf      = numel(f);

% Get number of spatial points 
nz      = numel(z);

% Define output frequency grid 3DM
f       = reshape(f(:),[1 1 nf]);

% z1 row vector
% eps is added in order to avoid the singularity of the coh
% function when z = 0
z1      = z(:)';
% Get number of structural points
nz      = numel(z);

% z2 column vector
z2      = z1';

% Define spatial grid 3DM
zz1     = z1(ones(nz,1),:,ones(nf,1));
zz2     = z2(:,ones(nz,1),ones(nf,1));

Su    = spectra(zz1,zz2,repmat(f,[nz,nz,1]));


function Sf1f1 = int2sf(z1,psi,cpsdf)

% z1 we have to ensure that z1 is a row vector
z1   = z1(:)';

% Column vector
z2   = z1';

nz   = numel(psi);
nf   = size(cpsdf,3);

% eigenmode i row vector
psi_iz1 = psi(:,1)';

% eigenmode j column vector
psi_jz2 = psi(:,1);

% First method (double nested trapz)fastest
% Define spatial grid 3DM
% These statements are equivalent to reshape but faster
psi_izz1     = psi_iz1(ones(nz,1),:,ones(nf,1));
psi_jzz2     = psi_jz2(:,ones(nz,1),ones(nf,1));
integrand    = psi_izz1.*cpsdf.*psi_jzz2;
Iz1          = trapz(z2,integrand, 2);
Sf1f1        = trapz(z1,Iz1,1);

% reshape the output
Sf1f1        = squeeze(Sf1f1);





