%***************************************************
%* Implemented by Isabelle Nunes,                  *
%* supervisioned by Prof. Dr. Victor L. Yoshimura  *
%* february - december, 2021                       *
%***************************************************
%% INPUT -------------------------------------------
%        |  A       | dim(n,n)  | 
%        |  P       | dim(n,n)  |  
%        |  alpha   | dim(1,1)  | 
%        |  r       | dim(1,1)  |
%        |  theta   | dim(1,1)  |
%% -------------------------------------------------    
%% OUTPUT ------------------------------------------
%        | R_1 | dim(n,n)   |
%        | R_2 | dim(2n,2n) |
%        | R_3 | dim(2n,2n) |
% --------------------------------------------------


function [R_1, R_2, R_3] = stableRegion(A, P, alpha, r, theta) %stableRegion(theta, r, P, A, alpha)

    if nargin < 3
        error("D-stable: Expecting at least 3 arguments (A, P, alpha, r, theta)")
        return
    end  

    theta = theta * (2*pi/360);             % theta em radianos 

    Z = zeros(size(A,1), size(A,2));     
    R_1 = A*P + P*A' + 2*alpha*P;           % eq 20 - Chilali.
    R_2 = [-r*P, A*P; P*A' , -r*P];         % eq 21 - Chilali
    
    % eq 22 - Chilali
    R_3 = [sin(theta)*(A*P + P*A') , cos(theta)*(A*P - P*A') ; cos(theta)*(P*A' - A*P) , sin(theta)*(A*P + P*A') ]; 
end
