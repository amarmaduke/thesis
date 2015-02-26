function [E, Ebt, Est, Evst] = simpEnergy( v, state )
	beta = state.beta;
	delta = state.delta;
	len = state.len;
	gamma = state.gamma;
	sigma = state.sigma;
	epsi = state.epsi;
	
	[a, ~] = size(v);
	
	E = zeros(1,a);
	Ebt = zeros(1,a);
	Est = zeros(1,a);
	Evst = zeros(1,a);

	for t = 1:a
		len_dr_1 = sqrt( (v(t,1) - delta)^2 + (v(t,3) - 0)^2 );
		len_dr_2 = sqrt( (v(t,2) - v(t,1))^2 + (v(t,4) - v(t,3))^2 );
		dot_prdct = (v(t,2) - v(t,1))*(v(t,1) - delta) + (v(t,4) - v(t,3))*v(t,3);
	
		
		%Ebt(t) = beta*(len_dr_1*len_dr_2/dot_prdct - 1 + len_dr_1/v(t,3) - 1);
		
		
		%Est(t) = gamma*((len_dr_1-len)^2 + (len_dr_2-len)^2);
		
		
		Evst(t) = pi*sigma*epsi*(.2*(sigma/v(t,3))^10 - (sigma/v(t,3))^4 ...
				+ .2*(sigma/v(t,4))^10 - (sigma/v(t,4))^4 );
			
		E(t) = Ebt(t) + Est(t) + Evst(t);
	end


end

