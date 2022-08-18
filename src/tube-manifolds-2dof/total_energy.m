function e = total_energy(orbit, parameters)

%   total_energy computes the total energy of an input orbit (represented
%   as M x N with M time steps and N = 4, dimension of phase space for the
%   model) for the 2 DoF solute-solvent model with LJ term in the potential.
% 
%   Orbit can be different initial conditions for the periodic orbit of
%   different energy. When trajectory is input, the output energy is mean.
%

    coeff_pe = parameters(3:end);
    
    r1 = orbit(:,1);
    r2 = orbit(:,2);
    pr1 = orbit(:,3);
    pr2 = orbit(:,4);
    
    e = (pr1.^2/(2*parameters(1))) + (pr2.^2/(2*parameters(2))) ...
        + compute_v1(coeff_pe, r1) ...
        + coeff_pe(6)*(coeff_pe(7) - r2).^2 ...
        + coeff_pe(8)/((r2 - r1).^12); 
    
    if length(e) > 1 
        e = mean(e);
    end
        
    
end
