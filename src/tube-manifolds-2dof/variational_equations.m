function PHIdot = variational_equations(t, PHI, parameters)

%        PHIdot = variational_equations(t, PHI, parameters) ;
%
% This here is a preliminary state transition, PHI(t,t0),
% matrix equation based on...
%
%        d PHI(t, t0)
%        ------------ =  Df(t) * PHI(t, t0)
%             dt


    x(1:4) = PHI(17:20);
    phi  = reshape(PHI(1:16),4,4);

%     % The first order derivative of the Hamiltonian.
%     dVdx = - 2*par(4)*par(5)*exp(-par(5)*x(1))*(exp(-par(5)*x(1)) - 1) - ...
%         4*par(6)*par(5)*x(2)^2*(x(2)^2 - 1)*exp(-par(6)*par(5)*x(1)) ;
%     dVdy = 8*x(2)*(2*x(2)^2 - 1)*exp(-par(6)*par(5)*x(1));
%     
%     % The following is the Jacobian matrix 
%     d2Vdx2 = - ( 2*par(4)*par(5)^2*( exp(-par(5)*x(1)) - ...
%                 2.0*exp(-2*par(5)*x(1)) ) - ...
%                 4*(par(6)*par(5))^2*x(1)^2*(x(2)^2 - 1)*exp(-par(6)*par(5)*x(1)) );
%             
%     d2Vdy2 = 8*(6*x(2)^2 - 1)*exp( -par(5)*par(6)*x(1) );
%     
%     d2Vdydx = -8*x(2)*par(5)*par(6)*exp( -par(5)*par(6)*x(1) )*(2*x(2)^2 - 1);
%     
%     d2Vdxdy = d2Vdydx;    
% 
%     Df    =[  0,     0,    1/par(1),    0;
%               0,     0,    0,    1/par(2);
%               -d2Vdx2,  -d2Vdydx,   0,    0;
%               -d2Vdxdy, -d2Vdy2,    0,    0];


    Df = jacobian(x, parameters);
    
    phidot = Df * phi; % variational equation

    PHIdot        = zeros(20,1);
    PHIdot(1:16)  = reshape(phidot,16,1); 
    PHIdot(17:20) = solutesolventLJ2dof(t, x, parameters);
    
%     PHIdot(17)    = x(3)/par(1);
%     PHIdot(18)    = x(4)/par(2);
%     PHIdot(19)    = -dVdx; 
%     PHIdot(20)    = -dVdy;

end