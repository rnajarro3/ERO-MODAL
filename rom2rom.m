function romOut = rom2rom(rom,site,modelParameters)
% Computes coefficients A_ij, A_ijk

% Compute average wind velocity
U_j           = modelParameters.U_1;

% Ensure a column vector
U_j           = U_j(:);

% Ensure a column vector
D_j               = rom.D(:);
nz                = length(D_j);

% Compute sdof-site data
alpha             = pi/2 - rom.thetaG;
cd_j              = zeros(nz,1);

for i = 1:nz
    fcd                = rom.airfoil{i}.cd;
    cd_j(i)            = fcd(alpha(i));
end

% Compute integrands
f_3               = 0.5*site.density.*cd_j.*D_j;
f_3               = f_3(:,ones(1,1));
f_2               =     site.density.*cd_j.*D_j.*U_j;
f_2               = f_2(:,ones(1,1));

% Ensure a column vetor
r            = rom.r(:);
Psi_ij       = rom.r(:)/rom.r(end);
% Psi_ij       = rom.Psi_ij(:);
% Second order tensor f2
a_2        = permute(f_2.*Psi_ij,[2,3,1]);
b_2        = permute(Psi_ij,[3,2,1]);
a_ijk      = multiprod(a_2,b_2, [1,2]);
A_ij       = trapz(r,a_ijk,3);

% % Third order tensor f3
a_3        = permute(f_3.*Psi_ij,[2,3,4,1]);
b_3        = permute(Psi_ij,[3,2,4,1]);
c_3        = permute(Psi_ij,[3,4,2,1]);       
a_ijkl     = multiprod(b_3,c_3,[2,3]);
a_ijkl     = multiprod(a_3,a_ijkl,[1,2]);
A_ijk      = trapz(r,a_ijkl,4);


%% Define output
romOut = rom;
romOut.U_1      = U_j;
romOut.A_ij     = A_ij;
romOut.A_ijk    = A_ijk;



end