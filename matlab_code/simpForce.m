function out = simpForce( t, v, state )
	persistent fcount;
	persistent ltime;
	persistent ptime;
	
	
	if t == 0
		ltime = cputime;
		fcount = 0;
		ptime = 0;
	end
	fcount = fcount + 1;
	
	%% Timing
	if mod(fcount,1000) == 0
		ttime = cputime;
		elpsd = ttime - ltime;
		ptime = ptime + elpsd;
		pcnt = 100 * (t / state.tend);
		cmpltnTime = ( ( ptime * 100 / pcnt ) - ptime ) / 60;
		fprintf('%.2f %%, %.5f s, %.5f m; %.5f m \n',pcnt,elpsd,ptime/60,cmpltnTime);
		ltime = cputime;
	end
	
	
	beta = state.beta;
	delta = state.delta;
	len = state.len;
	gamma = state.gamma;
	sigma = state.sigma;
	epsi = state.epsi;

	len_dr_1 = sqrt( (v(1) - delta)^2 + (v(3) - 0)^2 );
	len_dr_2 = sqrt( (v(2) - v(1))^2 + (v(4) - v(3))^2 );
	dot_prdct = (v(2) - v(1))*(v(1) - delta) + (v(4) - v(3))*v(3);
	
	out = zeros(length(v),1);
	
	b1 = 0; b2 = 0; b3 = 0; b4 = 0;
	g1 = 0; g2 = 0;
	
	%% x_1
	%b1 = (v(1) - delta)*len_dr_2/(len_dr_1*dot_prdct);
	%b2 = -(v(2) - v(1))*len_dr_1/(len_dr_2*dot_prdct);
	%b3 = -(v(2)-2*v(1)+delta)*len_dr_1*len_dr_2/dot_prdct^2;
	%b4 = (v(1) - delta)/(len_dr_1*v(3));
	%g1 = 2*(v(1) - delta)*(len_dr_1-len)/len_dr_1;
	%g2 = -2*(v(2) - v(1))*(len_dr_2-len)/len_dr_2;
	
	%first = (g1 + g2);
	out(1) = out(1) - beta*(b1 + b2 + b3 + b4) - gamma*(g1 + g2);
	%% x_2
	%b1 = (v(2) - v(1))*len_dr_1/(len_dr_2*dot_prdct);
	%b2 = -(v(1) - delta)*len_dr_1*len_dr_2/dot_prdct^2;
	%g1 = 2*(v(2) - v(1))*(len_dr_2-len)/len_dr_2;
	%second = (g1);
	out(2) = out(2) - beta*(b1 + b2) - gamma*g1;
	%% y_1
	%b1 = v(3)*len_dr_2/(len_dr_1*dot_prdct);
	%b2 = -(v(4) - v(3))*len_dr_1/(len_dr_2*dot_prdct);
	%b3 = -(v(4) - 2*v(3))*len_dr_1*len_dr_2/dot_prdct^2;
	%b4 = (1/len_dr_1) - len_dr_1/v(3)^2;
	%g1 = 2*v(3)*(len_dr_1-len)/len_dr_1;
	%g2 = -2*(v(4) - v(3))*(len_dr_2-len)/len_dr_2;
	v1 = -pi*epsi*( 2*(sigma/v(3))^11 - 4*(sigma/v(3))^5 );
	third = (v1);
	out(3) = out(3) - beta*(b1 + b2 + b3 + b4) - gamma*(g1 + g2) - v1;
	%% y_2
	%b1 = (v(4) - v(3))*len_dr_1/(len_dr_2*dot_prdct);
	%b2 = -v(3)*len_dr_1*len_dr_2/dot_prdct^2;
	%g1 = 2*(v(4) - v(3))*(len_dr_2-len)/len_dr_2;
	v1 = -pi*epsi*( 2*(sigma/v(4))^11 - 4*(sigma/v(4))^5 );
	fourth = (v1);
	out(4) = out(4) - beta*(b1 + b2) - gamma*g1 - v1;
	%[first, second, third, fourth]
	%[third, fourth]
	%pause
end

