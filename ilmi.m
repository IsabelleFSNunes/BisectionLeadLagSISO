%***************************************************
%* Implemented by Isabelle Nunes,                  *
%* supervisioned by Prof. Dr. Victor L. Yoshimura  *
%* february - december, 2021                       *
%***************************************************
%% INPUT -------------------------------------------
%        |  A  | dim(n,n) | matrix with the eigenvalues of system
%        |  B  | dim(m,n) | matrix with the inputs of system
%        |  C  | dim(p,n) | matrix with the output of system
%        | rho | dim(1,1) | stop condition
%% -------------------------------------------------    
%% OUTPUT ------------------------------------------
%        | isFeasible | dim(1,1) | boolean
%        |  F         | dim(m,p) | matrix of stabilizing static output gain
%% --------------------------------------------------

function [isFeasible, F, gama, alpha2] = ilmi(ac, A, B, C, theta, rho)

    if nargin ~= 6 
        if nargin == 5 
            rho = 10E-3                 % value rho - default 
        else 
            error("ilmi: Expecting 5 or 6 arguments (ac, A, B, C, theta, rho), the variable rho is optional.")
            return
        end     % end if 02
    end     % end if 01

    % To verify the matrices compatibility
    m = size(B,2);
    n = size(A,1);
    p = size(C,1);

    if size(A,2) ~= n || size(B,1) ~= n || size(C,2) ~= n
        error('ilmi: The Matrices are incompatible.')
        return
    end

    % adjusting the settings        
    op  = sdpsettings('solver', 'sdpt3');
    op2 = sdpsettings('solver','bisection','bisection.solver', 'sdpt3', 'sdpt3.rmdepconstr', 1);

    % requirement
    pkg load control 

    % initilize some local variables 
    index = 1;                              % initialize i, here called of index
    can_continue = true;                    % flag
    isFeasible = false;                     % flag for function return 
    alpha = sdpvar(1);
    F = sdpvar(size(B,2), size(C,1));       % initialize F
    Q = eye(size(A,1), size(A,2));          % selection of Q 

    %% Starts Algorithm ILMI - Cao

    % Step 01.
    % apply Equation's Riccati to find first X
    X = care(A, B, Q, eye(size(B,2)) );                          
    P = sdpvar(size(A,1), size(A,2) );              

    while can_continue
        printf("ilmi: iteration i = %d\n ", index);

        % Step 02. 
        % elements of M -------------------------------------------
        m_11 = A'*P + P*A - X*B*B'*P - P*B*B'*X + X*B*B'*X - alpha*P;   
        m_12 = (B'*P + F*C)';

        m_21 = B'*P + F*C;
        m_22 = -eye(size(B, 2), size(B, 2));

        M = [m_11 m_12; m_21 m_22];
        % ---------------------------------------------------------

        %  norm H infinity
        gama = sdpvar(1);              % gama as variable
        H_inf = h_inf(P, F, gama, A, B, C);
        
        %  D-stable region
        alpha2 = sdpvar(1);            % alpha eq 21 Chilali
        ratio = 50*abs( min(eig(A)) ); % defining a limitant 
        [R_1, R_2, R_3] =  stableRegion(A+B*F*C, P, alpha2, r, theta);

        %  subject to ... constraints list
        %  mandatory:    M <= 0, P >= 0 
        %  optional :    H_inf <= 0, R_1 <= 0, R_2 <= 0, R_3 <= 0
        LMI = [ M <= 0, P >= 0, H_inf <= 0, R_1 <= 0];
        
        diagnosticBisection = bisection( LMI, alpha, op2); 
        
        printf("Gama = %f \t alpha2 = %f\n", value(gama), value(alpha2) );
        value(F)
        if diagnosticBisection.problem == 0
            printf("\n\nThis LMI is feasible with alpha: %f\n", value(alpha));

            % Step 03. To verify if the system is stable with this alpha 
            if value(alpha) <= 0
                printf("\t F is stabilizing static output gain\n");
                can_continue = false;
                isFeasible = true;
            else 
                printf("\t alpha > 0 --> F is NOT stabilizing gain. \n");

                % Step 04. Minimize trace(P) 
                F = sdpvar(size(B,2), size(C,1) );   % rebooting F

                % update the elements of M, using P and alpha found in bisection
                m_11 = A'*P + P*A - X*B*B'*P - P*B*B'*X + X*B*B'*X - value(alpha)*P;
                m_21 = B'*P + F*C;
                m_22 = -eye(size(B, 2), size(B, 2));
                
                M = [m_11, m_21'; m_21, m_22];

                LMI = [M <= 0, P >= 0];      % constraints list 
                
                diagnosticFeasible = optimize( LMI, trace(P), op);       

                % 0 = Successfully solved
                if diagnosticFeasible.problem == 0                     

                    % Step 05. 
                    if norm( X - value(P) ) >= rho
                        % Prepare the switch next iteration
                        % update X
                        X = value(P);                      
                        
                        % rebooting F
                        F = sdpvar(size(B,2), size(C,1) ); 
                        
                        % rebooting alpha
                        alpha = sdpvar(1);  
                        
                        % update iteration
                        index++; 
                        
                        % go to step 2.
                        can_continue = true;               

                    % Step 06. 
                    else    
                        printf("Infeasible - Step 06\n");
                        can_continue = false;
                    end     % end if and else 
                % 1 - Infeasible problem
                elseif diagnosticFeasible.problem == 1                
                    printf("Step 04.2 ");
                    can_continue = false;
                else
                    printf("Infeasible - Step 04\n");
                    can_continue = false;
                end     % end if and else

            end     % end if (and else) - 02
        
        else
            printf("Infeasible Step 02\n");
            can_continue = false;
        end     % end if (and else)  - 01

    end     % end while 
    
    alpha = value(alpha);
end     % end function 