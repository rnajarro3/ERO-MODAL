 function  dydt = odeRomModel(t,y,rom)
 
dydt            = zeros(size(y));

% Time computation of external force (ft should be a column vector)
fj              = rom.fj;
tj              = rom.t(:);
ft              = interp1(tj,fj,t);
ft              = ft(:);

% Time computation of external force due to fluctuation times to structural
% velocity
Ft_ij           = rom.Ft_jk;
F_ij            = squeeze(interp1(tj,Ft_ij,t));

% slice y-vector according to physical interpretation based on number of
% modes
z_1             = y(1);
z_2             = y(2);

% Get rom entities
omega2         = rom.wT^2;
cs             = 2 * rom.zeta * rom.wT;
Inv_m_ij       = rom.r(end)/rom.IT;
% Inv_m_ij        = rom.r(end)^2/rom.IT;
A_ij           = rom.A_ij;
A_ijk          = rom.A_ijk;

%% Tensorial products

r_1              = z_1';
A_ijk_zjzk       = r_1*reshape(...
                   r_1*reshape(...
                   shiftdim(A_ijk,1),...
                   1, 1),...
                   1, 1);
A_ijk_zjzk       = A_ijk_zjzk';

% Structural dynamic ode
dydt(1)      = Inv_m_ij*(ft - F_ij*z_1 - A_ij*z_1 + A_ijk_zjzk) ...
                  - cs*z_1 ...
                  - omega2*z_2;

% Relation between velocities and displacements
dydt(2) = z_1;

end
