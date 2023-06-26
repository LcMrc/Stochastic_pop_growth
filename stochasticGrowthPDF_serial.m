function PN = stochasticGrowthPDF_serial(N0, K, b, g, t_th, GrowthModel)

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

    PN = NaN(K+1, length(t_th));
    PN(1 : N0, :) = 0;

	nd = 128; % number of digits in the rounding off of rates (prevents unique from distinguishing rates different by double precision only)
	digits(nd)

    L = [];
    
    for N = N0 : K-1
        
    	iN = N+1;

        L = [L vpa(birthRate(N))];

        n = length(L); % total number of rates
		l = unique(L,'stable'); % get list of independent rates
		n_rate = length(l); % number of unique rates

		P = NaN(1, length(t_th));

		switch n_rate
            
			% Case 1: n identical rates (=single rate)
			case 1
                
				P = l^(n-1)/factorial(n-1)*t_th.^(n-1).*exp(-l*t_th);
				mc = length(L);

			% Case 2: n distinct rates
			case n 
                
				A1 = repmat(l',[1,n]); A2 = repmat(l',[1 n]) - repmat(l,[n 1]); A2 = A2 + diag(l) - diag(A2); A = A1./A2;

                for iT = 1 : length(t_th)
                    
					tt = t_th(iT);
					P(iT) = sum(l.*exp(-tt*l).*prod(A));

                end
                
				P = P / L(end);
				mc = 1;
			
			% Case 3: n rates, some equal, others distinct
            otherwise 
                
				% Get degeneracy of each rate
				c = cell2mat(arrayfun(@(x)length(find(L == x)), l, 'Uniform', false));
				mc = max(c);
                
				if(mc > 2)
                    
					error('Highest rate degeneracy is larger than 2, not yet implemented!');
                    
                end
                
				A = 0; 
                
				for ii = 1 : n_rate
					
                    ci = c(ii);
                    
					for kk = 0 : ci-1
                        
						A = A + MPDCoeffs(ii,kk,l,c)/(factorial(kk)*factorial(ci-kk-1))*t_th.^(ci-kk-1).*exp(-l(ii)*t_th);
					
                    end
                    
                end
                
				P = prod(l.^c)*A;
				P = P / L(end);
                
		end

		if( max(abs(P)) > 1 )
            
			error(['Probability exceeded 1 for N = ' num2str(N) ' with max(P) = ' num2str(max(abs(P)))])
            
		end

		PN(iN,:) = P;

		disp(['Done with N = ' num2str(N) ', found n_rate = ' num2str(n_rate) ' out of ' num2str(length(L)) ' with max deg. = ' num2str(mc)]);
    
    end

    % Use normalization of probability for the case N=K
    PN(K+1, :) = 1-sum(PN(1 : K,:));
    
end

function MPD = MPDCoeffs(ii,kk,l,c)

	ll = [l(1:ii-1) l(ii+1:end)]; 
	cc = [c(1:ii-1) c(ii+1:end)]; 

	switch kk
        
		case 0
			MPD = prod((ll - l(ii)).^(-cc)); 
		
		case 1
			MPD = 0; 
            
			for jj = 1:length(ll)
                
				MPD = MPD + (-cc(jj))*prod((ll - l(ii)).^(-cc))/(ll(jj) - l(ii));
                
			end

    end

end