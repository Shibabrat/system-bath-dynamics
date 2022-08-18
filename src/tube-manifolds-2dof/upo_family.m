function [x0po,T] = upo_family(eqNum,Ax1,Ax2,nFam,po_fam_file,parameters) 

% Generate a family of periodic orbits (po) given a pair of seed initial
% conditions and periods
%
% Shibabrat Naik (modified: 12-Mar-20)

    % delt = guessed change in period between successive orbits in family

    delt = -1.e-9 ;   % <==== may need to be changed
%     delt = -1.e-12 ;   

    N = 4 ; % dimension of phase space

    x0po = zeros(nFam,N) ;
    T    = zeros(nFam,1) ;
    energyPO = zeros(nFam,1) ;
    
    [x0poGuess1,TGuess1] = upo_guess_linear(Ax1,parameters);
    [x0poGuess2,TGuess2] = upo_guess_linear(Ax2,parameters);

    % Get the first two periodic orbit initial conditions
    iFam = 1 ;
    FAMNUM = sprintf('::poFamGet : number %d',iFam) ;
    disp(FAMNUM);
    [x0po1,tfpo1] = upo_diff_corr_fam(x0poGuess1', parameters);  % see below 
    energyPO(iFam) = total_energy(x0po1, parameters) ;


    iFam = 2 ;
    FAMNUM = sprintf('::poFamGet : number %d',iFam) ;
    disp(FAMNUM) ;
    [x0po2,tfpo2] = upo_diff_corr_fam(x0poGuess2', parameters) ;               
    energyPO(iFam) = total_energy(x0po2, parameters) ;


    x0po(1:2,1:N) = [x0po1(:)'  ; x0po2(:)' ];
    T   (1:2)     = [2*tfpo1    ; 2*tfpo2     ]; 

    %Generate the other members of the family using numerical continuation
    for iFam = 3:nFam

        FAMNUM = sprintf('::poFamGet : number %d', iFam) ;
        disp(FAMNUM) ;
        
        dx  = x0po(iFam-1,1) - x0po(iFam-2,1) ;
        dy  = x0po(iFam-1,2) - x0po(iFam-2,2) ;
%         dOrbit = x0po(iFam-1,:) - x0po(iFam-2,:);
        dt  = T(iFam-1) - T(iFam-2) ;

%         x0po_g = x0po(iFam-1,:)' + dOrbit';
        t1po_g =   (T(iFam-1) + dt)/2 + delt ;
        x0po_g = [ (x0po(iFam-1,1) + dx) (x0po(iFam-1,2) + dy) 0 0] ;
%         t1po_g = T(iFam-2) + abs(dt);

      % differential correction takes place in the following function
        [x0po_iFam,tfpo_iFam] = upo_diff_corr_fam(x0po_g, parameters) ; 

        x0po(iFam,1:N) = x0po_iFam ;
%         T   (iFam)     = 2*t1_iFam ;	
        T(iFam)        = 2*tfpo_iFam;

        energyPO(iFam) = total_energy(x0po(iFam,:), parameters) ;
        
        if mod(iFam,10) == 0
            dum = [x0po T] ;
    %     save x0po_T.dat -ascii -double dum
        end

    end

    dum = [x0po T energyPO] ;
    save(po_fam_file,'dum','-ascii','-double');

end















