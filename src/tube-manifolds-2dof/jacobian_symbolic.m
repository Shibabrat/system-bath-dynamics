function Df = jacobian_symbolic(eqPt, parameters)

%%% par = [MASS_A MASS_B EPSILON_S D_X LAMBDA ALPHA];

%     global OMEGA_X OMEGA_Y DELTA
    
    syms r1 r2 pr1 pr2 
    
    coeff_pe = parameters(3:end);
    xe = eqPt(1);
    ye = eqPt(2);
    vxe = eqPt(3);
    vye = eqPt(4);
    
    dv1_dr1  = (2-1)*coeff_pe(2)*(r1.^(2-2)) + (3-1)*coeff_pe(3)*(r1.^(3-2)) ...
        + (4-1)*coeff_pe(4)*(r1.^(4-2)) + (5-1)*coeff_pe(5)*(r1.^(5-2));
    dvint_dr1 = 12*coeff_pe(8)/((r2 - r1).^13);
    
    dv2_dr2 = - 2*coeff_pe(6)*(coeff_pe(7) - r2);
    dvint_dr2 = -12*coeff_pe(8)/((r2 - r1).^13);
    
    %%% Use vector differentiation 
    f(r1, r2, pr1, pr2) = [pr1/(parameters(1)); 
                            pr2/(parameters(2)); 
                            - ( dv1_dr1 + dvint_dr1 ) ; 
                            - ( dv2_dr2 + dvint_dr2 ) ];
    
    DfMat(r1, r2, pr1, pr2) = jacobian(f,[r1 r2 pr1 pr2]);
    
    Df = double(DfMat(xe,ye,vxe,vye));
    
end