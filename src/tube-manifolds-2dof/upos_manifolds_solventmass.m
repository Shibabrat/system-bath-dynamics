% Script for continuation of the unstable periodic orbits in the mass of
% the solvent parameter. 
% Setting up parameters and global variables
global eqNum deltaE

MASS_SOLUTE = 1.0; 
% coeff_pe = [321.904484, -995.713452, 1118.689753, -537.856726, ...
%         92.976121, 1.0, 1.0, 0.00]; 
coeff_pe = [321.904484, -995.713452, 1118.689573, -537.856726, ...
            92.976121, 1.0, 1.0, 0.01];
        
for MASS_SOLVENT = 1.0 %:0.1:3.0

    parameters = [MASS_SOLUTE, MASS_SOLVENT, coeff_pe];

    eqNum = 2;  
    [eqPt] = equilibrium_pt(eqNum, parameters);

    % energy of the saddle equilibrium point
    eSaddle = total_energy(eqPt, parameters);

%     target_energy = 3.5;
    
    target_energy = 3.691966889;
%     target_energy = 2.1;

    %%% 

    nFam = 100; % use nFam = 10 for low energy

    % first two amplitudes for continuation procedure to get p.o. family
    % Ax1  = 2.e-5; % initial amplitude (1 of 2) values to use: 2.e-3
    % Ax2  = 2*Ax1; % initial amplitude (2 of 2)

    % amplitudes for mu_2 = 0.1 to 1.0
    Ax1 = 1.e-8;
    Ax2 = 2*Ax1;
    
    % amplitudes for uncoupled
%     Ax1 = 1.e-5;
%     Ax2 = 2*Ax1;

    tic;
    %  get the initial conditions and periods for a family of periodic orbits
    po_fam_file = ['x0_tp_fam_eqPt',num2str(eqNum), ...
                    'mass_solvent',num2str(num2str(MASS_SOLVENT)), ...
                    '_solutesolventLJ2dof.txt'];
    [po_x0Fam,po_tpFam] = upo_family(eqNum, Ax1, Ax2, nFam, po_fam_file, parameters) ; 

    poFamRuntime = toc;

    x0podata = [po_x0Fam, po_tpFam];


    deltaE = target_energy - eSaddle;

    po_fam_file = ['x0_tp_fam_eqPt',num2str(eqNum), ...
                    'mass_solvent',num2str(num2str(MASS_SOLVENT)), ...
                    '_solutesolventLJ2dof.txt'];
    %     target_energy = eSaddle + deltaE; 
    fprintf('Loading the periodic orbit family from data file %s \n',po_fam_file); 
    x0podata = importdata(po_fam_file);

    po_brac_file = ['x0po_T_energyPO_eqPt',num2str(eqNum), ...
                    'mass_solvent',num2str(num2str(MASS_SOLVENT)), ...
                    '_brac',num2str(deltaE),'_solutesolventLJ2dof.txt'];
    tic;

    % [x0poTarget,TTarget] = bracket_POEnergy_bp(target_energy, x0podata, po_brac_file);
    [x0poTarget,TTarget] = upo_bracket_target_energy(target_energy, ...
                            x0podata, po_brac_file, parameters);
    poTarE_runtime = toc;

    save(['model_parameters_eqPt',num2str(eqNum), ...
            'mass_solvent',num2str(num2str(MASS_SOLVENT)), ...
            '_E',num2str(deltaE), ...
            '_solutesolventLJ2dof.txt'], 'parameters', '-ASCII', '-double');


    % target specific periodic orbit
    % Target PO of specific energy with high precision; does not work for the
    % model 

    po_target_file = ['x0po_T_energyPO_eqPt',num2str(eqNum), ...
                        'mass_solvent',num2str(num2str(MASS_SOLVENT)), ...
                        '_DelE',num2str(deltaE),'_solutesolventLJ2dof.txt'];

    [x0_PO, T_PO, e_PO] = upo_bisect_target_energy(x0poTarget, ...
                            target_energy,po_target_file,parameters);


    data_path = ['./x0po_T_energyPO_eqPt', num2str(eqNum), ...
                'mass_solvent',num2str(num2str(MASS_SOLVENT)), ...
                '_DelE', num2str(deltaE), '_solutesolventLJ2dof.txt'];


end


%%

