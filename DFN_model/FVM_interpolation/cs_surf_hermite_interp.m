function cs_surf = cs_surf_hermite_interp(cs,r_index_Nr,cs_outer_ghost,param,domain)
%   cs_surf_hermite_interp computes the surface particle concentration
%   using 3rd order hermite interpolation.
%
%   Inputs:
%       cs: noramlized concentration along particle radius direction [-]
%       r_index_Nr: r index at r = Nrp
%       cs_outer_ghost: normalized particle concentration at surface outer ghost node [-]
%       param: Parameter structure
%       domain: string letter indicating electrode domain
%           -> 'p': positive electrode
%           -> 'n': negative electrode
%
%   Output:
%       cs_surf: noramlized surface particle concentration [-]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% COBRAPRO: Co-simulation Battery Modeling for Accelerated Parameter Optimization

% Copyright (c) 2024 CO-simulation BatteRy modeling for Accelerated PaRameter Optimization (COBRAPRO)
% COBRAPRO is freely distributed under the MIT License 
% v1.0.0.: Released March, 2024 

% Main contributors: 
% Sara Ha (sungyeon.sara.ha@stanford.edu)
% Simona Onori (sonori@stanford.edu)
% Stanford Energy Control Group, Stanford University
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

switch domain
    case 'p'
        r = param.rp;
        Nr = param.Nrp;
        Rp = param.Rpp;
    case 'n'
        r = param.rn;
        Nr = param.Nrn;
        Rp = param.Rpn;
end

r_Nr = r*(Nr-0.5); % radial distance from r = 0 to r = Nrp or r = Nrn
r_Nr_1 = r*(Nr-1.5); % radial distance from r = 0 to r = Nrp-1 or r = Nrn-1
p1 = (1+2*(Rp-r_Nr_1)/(r_Nr-r_Nr_1))*cs(r_index_Nr-1)+(Rp-r_Nr_1)*((cs(r_index_Nr)-cs(r_index_Nr-2))/(2*r));
p2 = ((Rp-r_Nr)/(r_Nr_1-r_Nr))^2;
p3 = (1+2*(r_Nr-Rp)/(r_Nr-r_Nr_1))*cs(r_index_Nr)+(Rp-r_Nr)*((cs_outer_ghost-cs(r_index_Nr-1))/(2*r));
p4 = ((Rp-r_Nr_1)/(r_Nr-r_Nr_1))^2;
cs_surf = p1*p2 + p3*p4;
