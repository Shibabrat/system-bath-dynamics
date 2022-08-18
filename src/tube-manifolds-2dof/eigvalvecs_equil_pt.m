function [Es,Eu,Ec,Vs,Vu,Vc] = eigvalvecs_equil_pt(equil_pt, parameters)

%   [Es,Eu,Ec,Vs,Vu,Vc] = eigvalvecs_equil_pt(ep,mu);
%
%   Find all the eigenvectors locally spanning the 3 subspaces for the phase
%   space in an infinitesimal region around an equilibrium point 
%
%   Our convention is to give the +y directed column eigenvectors
%
%   NOTE: Uses Matlab's built-in numerical algorithms to solve the eigenvalue
%       equation which arises
%

    Df = jacobian(equil_pt, parameters)

    [Es,Eu,Ec,Vs,Vu,Vc] = eigGet(Df,0);	% find the eigenvectors
                                        % i.e., local manifold approximations 

    % give the +y directed column vectors 
%     if Vs(2)<0, Vs=-Vs; end
%     if Vu(2)<0, Vu=-Vu; end

end




