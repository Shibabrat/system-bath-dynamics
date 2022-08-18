mass_solvent = 0.1;
total_energy = 3.691966889;
alpha = 0.4;
fs = draw_energysurf(mass_solvent, total_energy, alpha);
% hold on
% mass_solvent = 10;
% fs = draw_energysurf(mass_solvent, total_energy, alpha);


%% Test integration 

tspan = [0, 2];
Ax1  = 1.e-4;

MASS_SOLUTE = 1.0; 
coeff_pe = [321.904484, -995.713452, 1118.689753, -537.856726, ...
        92.976121, 1.0, 1.0, 0.01]; 

MASS_SOLVENT = 0.1;     % RGM, BC, SW (2019)
parameters = [MASS_SOLUTE, MASS_SOLVENT, coeff_pe];

eqNum = 2;  
[eqPt] = equilibrium_pt(eqNum, parameters);

[x0poGuess1, TGuess1] = upo_guess_linear(Ax1,parameters);
% x0po = load('x0po_T_energyPO_eqPt2_DelE0.21872_solutesolventLJ2dof.txt');

RelTol = 3.e-14; AbsTol = 1.e-14;
% options = odeset('RelTol',RelTol,'AbsTol',AbsTol,'Events', @half_period_event);
options = odeset('RelTol',RelTol,'AbsTol',AbsTol);

[tt,xx] = ode45(@(t, y)solutesolventLJ2dof(t, y, parameters), tspan, ...
                x0poGuess1, options);

critical_energy = total_energy(eqPt, parameters);

%%

numpts = 200;
xVec = linspace(0.7, 2.3, numpts);
yVec = linspace(1.4, 5.15, numpts);
[xMesh, yMesh] = meshgrid(xVec, yVec);
peMesh = potential_energy(xMesh, yMesh, parameters);

pe_range = -7:0.1:6;
target_energy = 3.691966889;
contour(xMesh, yMesh, peMesh, pe_range)
hold on
contour(xMesh, yMesh, peMesh, [target_energy target_energy], '-k')
colormap winter
colorbar
hold on
% plot(xx(:,1), xx(:,2), '-r','Linewidth', 2)
% plot(xx(1,1), xx(1,2), 'xk')

%%




