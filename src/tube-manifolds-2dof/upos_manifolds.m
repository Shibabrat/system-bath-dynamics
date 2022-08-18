global eqNum deltaE

% Setting up parameters and global variables
N = 4;          % dimension of phase space
MASS_SOLUTE = 1.0; 

% coeff_pe = [321.904484, -995.713452, 1118.689753, -537.856726, ...
%         92.976121, 1.0, 1.0, 0.01]; 
coeff_pe = [321.904484, -995.713452, 1118.689573, -537.856726, ...
            92.976121, 1.0, 1.0, 0.01];

MASS_SOLVENT = 1.0;

parameters = [MASS_SOLUTE, MASS_SOLVENT, coeff_pe];

eqNum = 2;  
[eqPt] = equilibrium_pt(eqNum, parameters);

% energy of the saddle equilibrium point
eSaddle = total_energy(eqPt, parameters);


%% 

nFam = 100; % use nFam = 10 for low energy

% first two amplitudes for continuation procedure to get p.o. family
% Ax1  = 2.e-5; % initial amplitude (1 of 2) values to use: 2.e-3
% Ax2  = 2*Ax1; % initial amplitude (2 of 2)

% amplitudes for mu_2 = 0.1
Ax1 = 1.e-8;
Ax2 = 2*Ax1;

% Ax1 = 1.e-5;
% Ax2 = 2*Ax1;

tic;
%  get the initial conditions and periods for a family of periodic orbits
po_fam_file = ['x0_tp_fam_eqPt',num2str(eqNum), ...
                '_solventmass',num2str(MASS_SOLVENT), ...
                '.txt'];
[po_x0Fam,po_tpFam] = upo_family(eqNum, Ax1, Ax2, nFam, po_fam_file, parameters) ; 

poFamRuntime = toc;

x0podata = [po_x0Fam, po_tpFam];


%%

target_energy = 3.691966889;
% paramstarttime = tic;    
    
deltaE = target_energy - eSaddle;

% po_fam_file = ['x0_tp_fam_eqPt',num2str(eqNum), ...
%                 '_solventmass',num2str(MASS_SOLVENT), ...
%                 '.txt'];
%     target_energy = eSaddle + deltaE; 
fprintf('Loading the periodic orbit family from data file %s \n',po_fam_file); 
x0podata = importdata(po_fam_file);

po_brac_file = ['x0po_T_energyPO_eqPt',num2str(eqNum), ...
                '_solventmass',num2str(MASS_SOLVENT), ...
                '_brac_deltaE',num2str(deltaE),'.txt'];
tic;

% [x0poTarget,TTarget] = bracket_POEnergy_bp(target_energy, x0podata, po_brac_file);
[x0poTarget,TTarget] = upo_bracket_target_energy(target_energy, ...
                        x0podata, po_brac_file, parameters);
poTarE_runtime = toc;

save(['model_parameters_eqPt',num2str(eqNum), ...
        '_solventmass',num2str(num2str(MASS_SOLVENT)), ...
        '_deltaE',num2str(deltaE), ...
        '.txt'], 'parameters', '-ASCII', '-double');




%% %

% target specific periodic orbit
% Target PO of specific energy with high precision; does not work for the
% model 

po_target_file = ['x0po_T_energyPO_eqPt',num2str(eqNum), ...
                    '_solventmass',num2str(MASS_SOLVENT), ...
                    '_deltaE',num2str(deltaE),'.txt'];

[x0_PO, T_PO, e_PO] = upo_bisect_target_energy(x0poTarget, ...
                        target_energy,po_target_file,parameters);


data_path = ['./x0po_T_energyPO_eqPt',num2str(eqNum), ...
            '_solventmass',num2str(MASS_SOLVENT), ...
            '_deltaE', num2str(deltaE), '.txt'];

%%

n_mfd_traj = 50;
frac = 0;
del = 1e-8;

x0po = importdata(data_path);

TPOFam = x0po(:,5); 
ePOFam = x0po(:,6);
nMed = size(x0po,1);
tmfd = 2.7; %2.7
% tmfd = 3*TPOFam(nMed);
% tmfd = 5.5*TPOFam(nMed);
% tmfd = 7.5*TPOFam(nMed);

% manistarttime = tic;    

stbl = 1;
dir = 1; % dir = 1: towards the top well with x as vertical axis
[xW,x0W] = upo_mani_globalize(x0po(nMed,1:4), TPOFam(nMed), frac, stbl, dir, ...
                            del, tmfd, n_mfd_traj, parameters);

% maniendtime = toc;

% maniruntime = maniendtime - manistarttime;

energyTube = ePOFam(nMed) ;
title(['Total energy: ', num2str(energyTube)]);

% end

paramendtime = toc;

% paramruntime = paramendtime - paramstarttime;


% % %

% manistarttime = tic;    

stbl = -1;
dir = 1; % dir = 1: towards the top well with x as vertical axis
[xW,x0W] = upo_mani_globalize(x0po(nMed,1:4),TPOFam(nMed), frac, stbl, dir, ...
                            del,tmfd,n_mfd_traj,parameters);

% maniendtime = toc;

% maniruntime = maniendtime - manistarttime;

energyTube = ePOFam(nMed) ;
title(['Total energy: ', num2str(energyTube)]);

% end


paramendtime = toc;


%%
hold on

numpts = 200;
xVec = linspace(0.7, 2.3, numpts);
yVec = linspace(1.45, 5.15, numpts);
[xMesh, yMesh] = meshgrid(xVec, yVec);
peMesh = potential_energy(xMesh, yMesh, parameters);

% target_energy = 3.691966889;
contour(yMesh, xMesh, peMesh, [target_energy target_energy], '-k', ...
        'linewidth', 2)
% view(90,-90)
set(gca,'FontSize', 40)
box on
grid off
xlabel('$y$', 'interpreter','latex', 'fontsize', 60)
ylabel(['$x' ...
    '$'], 'interpreter','latex', 'fontsize', 60)
title([])

% print -dpng -r300 stabletop_unstablebottom_solventmass_10e-1.png

%%



