function NN = GillespieStochasticGrowth(N0, K, b, g, Nit, t_sim, GrowthModel)

    switch GrowthModel
        
        case 'L' % Logistic
            
            birthRate = @(x) b*x*(1-x/K);
            
        case 'B' % Blumberg
            
            birthRate = @(x) b*x*(1-x/K)^g;
            
        case 'R' % Richards
            
            birthRate = @(x) b*x*(1-(x/K)^g);
            
        case 'G' % Gompertz
            
            birthRate = @(x) b*x*log(K/x);
            
        otherwise
            
            error('Growth Model not recognized!');
            
    end

    NN = NaN(length(t_sim), Nit);
    
    for i = 1 : Nit
        
        % Initialization
        N = N0;              % Initialization of the number of individuals
        t = 0;               % Initialization of time
        q = 1;
        NN(q, i) = N;
        q = q+1;

        while N ~= 0 && t < max(t_sim) && N ~= K
            
            % Compute the transition rates
            T = birthRate(N);   

            % Compute tau and then the new time
            r1 = rand;
            tau = 1/T*log(1/r1);
            t = t+tau;
            
            while t > t_sim(q) && q < length(t_sim)
                
                NN(q, i) = N;
                q = q+1;
                
            end

            N = N+1;
            
        end
    
        while q <= length(t_sim)
            
            NN(q, i) = N;
            q = q+1;
            
        end
        
    end
    
end