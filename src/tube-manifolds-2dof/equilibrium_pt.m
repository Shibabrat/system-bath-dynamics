function [eqPt] = equilibrium_pt(eqNum, parameters)

%equilibrium_pt solves the equilibrium points for the solute-solvent LJ 2
%DOF. These are numbered as going from left to right
%-------------------------------------------------------------------------------
%   Indices for the equilibrium points on the potential energy surface:
%
% 
%    center (EQNUM = 1)     saddle (EQNUM = 2)      center (EQNUM = 3)
%
%
%-------------------------------------------------------------------------------
%   
    
    % All the equilibrium points are solved numerically,
    if eqNum == 1 
        x0 = [1.0, 1.75, 0, 0];
    elseif 	eqNum == 2
        x0 = [1.5, 2.0, 0, 0];
%         x0 = [1.5, 0.5, 0, 0];  
    elseif 	eqNum == 3
        x0 = [2.0, 2.5, 0, 0];
    end
    
%%% F(xEq) = 0 at the equilibrium point, solve using in-built function
%    options = optimoptions('fsolve','Display','iter'); % Option to display output
%    [eqPt,fval] = fsolve(@func_vec_field_eq_pt, x0, options) % Call solver
    
    [eqPt, fval, ~] = ...
        fsolve(@(x) solutesolventLJ2dof(x, parameters), x0, ...
        optimset("jacobian","off")); % Call solver
 
end
function xDot = solutesolventLJ2dof(x, parameters)

    coeff_pe = parameters(3:end);
    xDot = zeros(length(x),1);
   
    
    xDot(1) = x(3)/parameters(1);
    xDot(2) = x(4)/parameters(2);

    dv1_dr1  = (2-1)*coeff_pe(2)*(x(1).^(2-2)) ...
        + (3-1)*coeff_pe(3)*(x(1).^(3-2)) ...
        + (4-1)*coeff_pe(4)*(x(1).^(4-2)) ...
        + (5-1)*coeff_pe(5)*(x(1).^(5-2));
    
    dvint_dr1 = 12*coeff_pe(8)/((x(2) - x(1)).^13);
    
    xDot(3) = - ( dv1_dr1 + dvint_dr1 );

    dv2_dr2 = - 2*coeff_pe(6)*(coeff_pe(7) - x(2));
    dvint_dr2 = -12*coeff_pe(8)/((x(2) - x(1)).^13);
    
    xDot(4) = - ( dv2_dr2 + dvint_dr2 );

end












