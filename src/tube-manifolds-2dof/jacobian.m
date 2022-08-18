function Df = jacobian(eqPt, parameters)
   
    
    coeff_pe = parameters(3:end);
    r1 = eqPt(:,1);
    r2 = eqPt(:,2);
    pr1 = eqPt(:,3);
    pr2 = eqPt(:,4);
    
    Df = zeros(4,4);

    
    
    d2v1 = (3 - 1)*(3 - 2)*coeff_pe(3)*(r1.^(3 - 3)) ...
            + (4 - 1)*(4 - 2)*coeff_pe(4)*(r1.^(4 - 3)) ...
            + (5 - 1)*(5 - 2)*coeff_pe(5)*(r1.^(5 - 3)) ;

    d2vint = (156*coeff_pe(end))/((r2 - r1).^14);
        
    deriv_wrt_r1_f3 = - d2v1 - d2vint;
    
    deriv_wrt_r2_f3 = d2vint;
    
    deriv_wrt_r1_f4 = d2vint;
    
%     deriv_wrt_r2_f4 = - 2*coeff_pe(end-3) - d2vint; % WTH!
    deriv_wrt_r2_f4 = - 2*coeff_pe(end-2) - d2vint; 

    Df(1,3) = 1/parameters(1);
    Df(2,4) = 1/parameters(2);
    
    Df(3,1) = deriv_wrt_r1_f3;

    Df(3,2) = deriv_wrt_r2_f3;
        
    Df(4,1) = deriv_wrt_r1_f4;
    
    Df(4,2) = deriv_wrt_r2_f4;

    
    

%     dv1_dr1  = (2-1)*coeff_pe(2)*(r1.^(2-2)) + (3-1)*coeff_pe(3)*(r1.^(3-2)) ...
%         + (4-1)*coeff_pe(4)*(r1.^(4-2)) + (5-1)*coeff_pe(5)*(r1.^(5-2));
%     dvint_dr1 = 12*coeff_pe(8)/((r2 - r1).^13);
%     
%     dv2_dr2 = - 2*coeff_pe(6)*(coeff_pe(7) - r2);
%     dvint_dr2 = -12*coeff_pe(8)/((r2 - r1).^13);
%     
%     %%% Use vector differentiation 
%     f(r1, r2, pr1, pr2) = [pr1/(parameters(1)); 
%                             pr2/(parameters(2)); 
%                             - ( dv1_dr1 + dvint_dr1 ) ; 
%                             - ( dv2_dr2 + dvint_dr2 ) ];
%     
%     DfMat(r1, r2, pr1, pr2) = jacobian(f,[r1 r2 pr1 pr2]);
%     
%     Df = double(DfMat(xe,ye,vxe,vye));
    
end




