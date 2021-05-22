function romdr  = getRomDynamicResponse(rom,site,t,wsssp,ic,modelParameters,varargin)
% getRomDynamicResponse gets the dynamic response of a ROM in time domain.
% 
% ROMDR  = getEbbDynamicResponse(ROM,SITE,TSPAN,WSSSP,IC,MP) gets dynamic
% response of a reduce order model defined by its ROM information, ROM; where
% SITE is a structure with the site information, TSPAN the time vector 
% that is simulated and WSSSP the time series of the external action. 
% The output structure ROMDR is defined as follows
%     t:        []
%     fj:       []
%     q:        []
%     v:        []
%     a:        []
%   eps:        []
%  deps:        []
% d2eps:        []
%   fzt:        []

[nu_i,nz,ny,nt,na]   = size(wsssp);
romdr         = cell(na,1);

%% ROM initial condition
y0          = zeros(2,1);

%% Compute ROM parameters
rom = rom2rom(rom,site,modelParameters); 

%% Forced response computation
for i = 1:na          
    
    [rom_a,f]    = sp2romAction(rom,site,wsssp(:,:,:,:,i));
    
    % Ensure t is a column vector
    rom_a.t      = t(:);

    % Resolution of ROM differential equation 
    romdr{i}          = romSolverODE(t,rom_a,y0);
    romdr{i}.fzt      = f;

end

end