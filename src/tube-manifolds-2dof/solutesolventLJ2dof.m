function xDot = solutesolventLJ2dof(t, x, parameters)

    coeff_pe = parameters(3:end);
    xDot = zeros(length(x),1);
   
    
    xDot(1) = x(3)/parameters(1);
    xDot(2) = x(4)/parameters(2);

    dv1_dr1  = (2-1)*coeff_pe(2)*(x(1).^(2-2)) ...
        + (3-1)*coeff_pe(3)*(x(1).^(3-2)) ...
        + (4-1)*coeff_pe(4)*(x(1).^(4-2)) ...
        + (5-1)*coeff_pe(5)*(x(1).^(5-2));
    
    dvint_dr1 = (12*coeff_pe(8))/((x(2) - x(1)).^13);
    
    xDot(3) = - ( dv1_dr1 + dvint_dr1 );

    dv2_dr2 = - 2*coeff_pe(6)*(coeff_pe(7) - x(2));
    dvint_dr2 = -(12*coeff_pe(8))/((x(2) - x(1)).^13);
    
    xDot(4) = - ( dv2_dr2 + dvint_dr2 );

end