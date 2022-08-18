function fs = draw_energysurf(mu2, H_val, alpha)


    % plot properties
    axesFontName = 'factory';
    % axesFontName = 'Times New Roman';
    axFont = 15;
    textFont = 15;
    labelFont = 20;
    lw = 2;    
    set(0,'Defaulttextinterpreter','latex', ...
        'DefaultAxesFontName', axesFontName, ...
        'DefaultTextFontName', axesFontName, ...
        'DefaultAxesFontSize',axFont, ...
        'DefaultTextFontSize',textFont, ...
        'Defaultuicontrolfontweight','normal', ...
        'Defaulttextfontweight','normal', ...
        'Defaultaxesfontweight','normal');

    % Parameters of the Hamiltonian
    MASS1 = 1; % mass of particle 1
    MASS2 = mu2; % mass of particle 2
    coeff_pe = [321.904484, -995.713452, 1118.689753, -537.856726, ...
        92.976121, 1.0, 1.0, 0.01];  
    
%     v1 = sum(coeff_pe(1:5).*r1.^(0:4));
    
    esurf = @(r1,r2,pr2) (pr2.^2)/(2*MASS2) ...
        + compute_v1(coeff_pe,r1) ...
        + coeff_pe(6)*(coeff_pe(7) - r2).^2 ...
        + coeff_pe(8)/((r2 - r1).^12) ...
        - H_val;
    
%     rgb_col = [51/255 153/255 51/255]; % green
    lightGrey1   = [0.85 0.85 0.85];
    lightGrey2   = [0.7 0.7 0.7];

    darkGrey1  = [0.4 0.4 0.4];
    darkGrey2  = [0.2 0.2 0.2];
    rgb_col = darkGrey1;
    
    xi = 0.7; xf = 2.3;
    yi = 1.4; yf = 5.15;
    pxi = -20; pxf = 20;

    fs = fimplicit3(esurf,[xi xf yi yf pxi pxf],...
        'EdgeColor','none','MeshDensity',100,'FaceAlpha',alpha,...
        'FaceColor',rgb_col);
  
    xlabel('$r_1$','FontSize',labelFont,'Interpreter','Latex');
    ylabel('$r_2$','FontSize',labelFont,'Interpreter','Latex');
    zlabel('$p_{r_2}$','FontSize',labelFont,'Interpreter','Latex');
    light;
    
    fs.FaceAlpha = alpha;
%     fs.AmbientStrength = 0.8;
    
end





