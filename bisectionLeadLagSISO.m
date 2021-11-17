%***************************************************
%* Implemented by Isabelle Nunes,                  *
%* supervisioned by Prof. Dr. Victor L. Yoshimura  *
%* february - december, 2021                       *
%***************************************************
%% INPUT -------------------------------------------
%        |  A       | dim(n,n) | 
%        |  B       | dim(m,n) |
%        |  C       | dim(p,n) |
%        | op_test  | dim(1,1) | number test
%        | delta    | dim(1,1) | stop condition bisection
%% -------------------------------------------------    
%% OUTPUT ------------------------------------------
%   Lead lag compensator coefficients 
%        | alpha | dim(1,1) |
%        | beta0 | dim(1,1) |
%        | beta1 | dim(1,1) |
% --------------------------------------------------

function [alpha, beta0, beta1, gama] = bisectionLeadLagSISO(A, B, C, op_test, delta)

    if nargin ~= 5
        if nargin == 4  
            delta = 0.1 ;                % value delta - default 
        else 
            error("ilmi: Expecting 3 or 4 arguments (A, B, C, delta), the variable delta is optional.")
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

    % Step 1. Initializing some variables --------------------
    % ac details, first estimative
    ac = - abs( ( min(eig(A)) - max(eig(A)) ) / 2 );    
    ac_min = -10* abs( min(eig(A)) );    % lower bound 
    ac_max = -10E-9;                     % upper bound (It tends to zero)

    % bounds bisections
    nFeasible = 0;             % iteration number of feasibles
    nFeasibleMAX = 5;          % iteration number of feasibles

    % angle (in degree) to consider D-stable region
    theta = 30;                
    theta_min = 30;        
    theta_max = 60;
    
    index = 0; 
    GoOn = true;               % bisection lead lag
    GoOn2 = false;             % bisection angle (overshoot)

    % to extract the important values
    local_file_out = sprintf("./Tests/Datasets/lmi-com-M-P-Hinf-R1/ex0%d.csv", op_test);
    filenameOutput = fopen(local_file_out, "w");
    fprintf(filenameOutput, " nFeasible, alpha, beta0, beta1, gama, mi,alpha_2\n");

    % matrices in closed loop (and space state extended)
    Af = [A, B; zeros(1,size(A,2)), ac]
    Bf = [B, zeros(size(B,1), 1); 0, 1]
    Cf = [C 0]
    % --------------------------------------------------

    % Step 2. to amplify the lower bound 
    while GoOn 
        printf("\nBisection Lead Lag : iteration %d\n", nFeasible);
        % using the ilmi function with ac specified 
        [isFeasible, F, gama, alpha_2] = ilmi(ac, Af,Bf,Cf, theta, 10E-5); 
             
        if isFeasible
            
            % bisection angle (overshoot) 
            while GoOn2
                theta = (theta_max + theta_min) / 2;

                % using the ilmi function to ac and theta specified 
                [isFeasible2, F, gama, alpha_2] = ilmi(ac, Af,Bf,Cf, theta, 10E-5);        
                printf("Bisection angle (Mp) = %d , iteration: %d \n", theta, index);

                % adjusting the limitants
                if isFeasible2
                    theta_max = theta;
                else
                    theta_min = theta;
                end

                if abs(theta_min - theta_max) <= 1.0 
                    GoOn2 = false;
                end % end if 
            end
            
            ac_max = ac;
            ac = 2*ac;
            printf("loop 1 - ac = %f \n", ac);
            value(F)
            nFeasible++;
        else                            
            if ac_max == 0 
                ac_min = ac;
            end
            GoOn = false;
        end

        % to load the coefficients
        if nFeasible <= nFeasibleMAX
            alpha = -ac_max;
            beta1 = value(F(1));                    
            beta0 = value(F(2)) + alpha*beta1;     


            % to extract the important values ---------------
            fprintf(filenameOutput,"%d, %f, %f, %f, %f, %f, %f\n", nFeasible, alpha, beta0, beta1, value(gama), value(gama)^2, value(alpha_2));
            % -----------------------------------------------

            % limit nFeasible
            if nFeasible == nFeasibleMAX
                GoOn = false;
                printf("This position is far to the left\n");
                fclose(filenameOutput);
                return 
            end

        end

    end

    GoOn = true;

    while GoOn
        % calculate the midpoint between the upper and lower bound
        ac = (ac_max + ac_min) / 2;         
        % using the ilmi function with ac and theta specified 
        [isFeasible, F,gama, alpha_2] = ilmi(ac, Af,Bf,Cf, theta, 10E-4);        
        printf("loop 2 - ac = %f \n", ac);
        value(F)
        
        index = 0;
        while GoOn2
            theta = (theta_max + theta_min) / 2;

            % using the ilmi function to ac and theta specified 
            [isFeasible2, F, gama, alpha_2] = ilmi(ac, Af,Bf,Cf, theta, 10E-5);        
            printf("Bisection angle (Mp) = %d , iteration: %d \n", theta, index);

            % adjusting the limitants
            if isFeasible2
                theta_max = theta;
            else
                theta_min = theta;
            end

            if abs(theta_min - theta_max) <= 1.0 
                GoOn2 = false;
            end % end if 
        end

        % to adjust the bounds
        printf("Is Feasible: ")
        isFeasible
        if isFeasible                   
            ac_max = ac;
        else
            ac_min = ac;
        end
        
        % stop condition to go out the loop 
        if abs( ac_min - ac_max ) <= delta      
            alpha = -ac_max;
            beta1 = value(F(1));                    
            beta0 = value(F(2)) + alpha*beta1;      

            % to extract the important values ---------------
            fprintf(filenameOutput,"%d, %f, %f, %f, %f, %f, %f\n", nFeasible, alpha, beta0, beta1, value(gama), value(gama)^2, value(alpha_2));
            % -----------------------------------------------
            GoOn = false;
            printf("Finalizing the Bisection Lead Lag\n");
            fclose(filenameOutput); 
            return 
        end 
    end         % end loop while 02

end         % end function 