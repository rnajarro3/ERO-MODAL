function [rom,fzt] = sp2romAction(rom,site,uspwss,varargin)
% Computes f_ij, Ft_jk

% Wind tilde time series
uspwss    = squeeze(uspwss)';
[nt,nz]   = size(uspwss);

% Ensure a column vector
r                 = rom.r(:);
Psi_ij            = rom.r(:)/rom.r(end);
% Psi_ij            = rom.Psi_ij(:);
% Ensure a row vector
c_j               = rom.D(:)';
cij               = c_j(ones(nt,1),:);

% Compute sdof-site data
alpha             = pi/2 - rom.thetaG;
cd_j              = zeros(1,nz);
for i = 1:nz
    fcd                = rom.airfoil{i}.cd;
    cd_j(i)            = fcd(alpha(i));
end
cd_ij                   = cd_j(ones(nt,1),:);

% External force distribution 
fzt               = 0.5*site.density.*cd_ij.*cij.*uspwss.^2.*sign(uspwss);

%% Compute modal external force due to utilde_i
a_2        = permute(fzt,[1,3,2]);
b_2        = permute(Psi_ij,[3,2,1]);
a_ijk      = multiprod(a_2,b_2, [1,2]);
f_ij       = trapz(r,a_ijk,3);

% Compute average wind velocity
U_j            = mean(uspwss);

% Ensure a row vector
U_j            = U_j(:)';
u              = uspwss - U_j(ones(nt,1),:);
fu             = site.density.*cd_ij.*cij.*u;
fu             = fu(:,:,ones(1,1,1));
a_3            = permute(fu,[1,4,3,2]);
b_3            = permute(Psi_ij,[3,2,4,1]);
c_3            = permute(Psi_ij,[3,4,2,1]);   

f_ijkl         = multiprod(b_3,c_3,[2,3]);
f_ijkl         = multiprod(a_3,f_ijkl,[1,2]);
Ft_jk          = trapz(r,f_ijkl,4);

% Overload rom output
rom.fj       = f_ij;
rom.Ft_jk    = Ft_jk;


end