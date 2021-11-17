%***************************************************
%* Implemented by Isabelle Nunes,                  *
%* supervisioned by Prof. Dr. Victor L. Yoshimura  *
%* february - december, 2021                       *
%***************************************************
%% INPUT -------------------------------------------
%        |  P    | dim(n,n) |  
%        |  F    | dim(p,n) |  
%        | gama  | dim(1,1) | 
%        |  A    | dim(n,n) | matrix A (with the eingvalues)
%        |  B_2  | dim(m,n) | matrix B 
%        |  C_2  | dim(p,n) | matrix C  
%% -------------------------------------------------    
%% OUTPUT ------------------------------------------
%        | H_inf | dim(n+2,n+2) |
% --------------------------------------------------

function H_inf = h_inf(P, F, gama, A, B_2, C_2) 

    if nargin ~= 6
        error("H_infinity: Expecting 6 arguments (P, F, gama, A, B_2, C_2)")
        return
    end     % end if 01
    B_1 = [1; zeros(size(A,1)-1, 1)];
    C_1 = [1, zeros(1, size(A,1)-1)];
    D_11 = 1;
    %D_12 = 0;
    D_12 = [0, 0];
    D_21 = 0;
    
    P_bar = [ [P, zeros(size(A,1),2)]; [zeros(2,size(A,1)), eye(2)] ];
    A_bar = [ [A, B_1, zeros(size(A,1),1)]; [zeros(1,size(A,1)), -gama*0.5 , 0]; [ C_1, D_11, -0.5*gama] ];
    B_bar = [ B_2; zeros(1, size(D_12,2)); D_12 ];
    C_bar = [ C_2, D_21, 0];
    
    H_inf = (P_bar * B_bar * F * C_bar) + (P_bar * B_bar * F * C_bar)' + A_bar' * P_bar + P_bar * A_bar ; 
end 