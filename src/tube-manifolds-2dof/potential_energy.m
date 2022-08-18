function pot_energy = potential_energy(r1, r2, parameters)
   
    
    coeff_pe = parameters(3:end);
    
    pot_energy = compute_v1(coeff_pe, r1) ...
        + coeff_pe(6)*(coeff_pe(7) - r2).^2 ...
        + coeff_pe(8)./((r2 - r1).^12);
                
    
end


