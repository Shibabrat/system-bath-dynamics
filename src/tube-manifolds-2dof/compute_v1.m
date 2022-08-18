function v1 = compute_v1(coeff_pe,r1)

    v1 = coeff_pe(1) + coeff_pe(2)*r1 + coeff_pe(3)*r1.^2 ...
        + coeff_pe(4)*r1.^3 + coeff_pe(5)*r1.^4;
    

end