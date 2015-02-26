function [E, Ebt, Est, Evst, Evdwt] = energy(y, state, imap )

	n = state.n;
	m = length(n);
	N = sum(n);
	beta = state.beta;
	delta = state.delta;
	gamma = state.gamma;
	sigma = state.sigma;
	density = state.density;
	len = state.len;
	epsi = state.epsi;
	
	[a, ~] = size(y);
	
	E = zeros(1,a);
	Ebt = zeros(1,a);
	Est = zeros(1,a);
	Evst = zeros(1,a);
	Evdwt = zeros(1,a);
	
	for t = 1:a
		Eb = 0;
		Es = 0;
		Evs = 0;
		Evdw = 0;
		
		%%{
		z = 1;
		for j = 1:m
			for i = 1:n(j)
				xI = z; % x index for vector (j, i)
				yI = xI + N; % y index
				
				%%{
				if i == 1
					len_1 = len;
					len_2 = sqrt((y(t,xI)-delta(j))^2 + y(t,yI)^2);
					dotprdct = y(t,yI)*len;
					Eb = Eb + len_1*len_2/dotprdct - 1;
				elseif i == 2
					len_1 = sqrt((y(t,xI-1)-delta(j))^2 + (y(t,yI-1)-0)^2); 
					len_2 = sqrt((y(t,xI)-y(t,xI-1))^2 + (y(t,yI)-y(t,yI-1))^2);
					dotprdct = (y(t,xI-1)-delta(j))*(y(t,xI)-y(t,xI-1)) + ...
								(y(t,yI)-y(t,yI-1))*(y(t,yI-1)-0);
					Eb = Eb + len_1*len_2/dotprdct - 1;
				else
					len_1 = sqrt((y(t,xI-1)-y(t,xI-2))^2 + (y(t,yI-1)-y(t,yI-2))^2); 
					len_2 = sqrt((y(t,xI)-y(t,xI-1))^2 + (y(t,yI)-y(t,yI-1))^2);
					dotprdct = (y(t,xI-1)-y(t,xI-2))*(y(t,xI)-y(t,xI-1)) + ...
								(y(t,yI)-y(t,yI-1))*(y(t,yI-1)-y(t,yI-2));
					Eb = Eb + len_1*len_2/dotprdct - 1;
				end
				%}
				
				%%{
				if i == 1
					dist = sqrt((y(t,xI)-delta(j))^2 + (y(t,yI))^2);
					Es = Es + (dist - len)^2;
				else
					dist = sqrt((y(t,xI)-y(t,xI-1))^2 + (y(t,yI)-y(t,yI-1))^2);
					Es = Es + (dist - len)^2;
				end
				%}
				
				Evs = Evs + LJsubstrate_Energy(y(t,yI),epsi,sigma,density);
				
				z = z + 1;
			end
		end
		
		Eb = beta*Eb;
		Es = gamma*Es;

		%%{
		for i = 1:N
			for k = i+1:N
				x1_I = i;
				y1_I = x1_I + N;
				x2_I = k;
				y2_I = x2_I + N;
				
				v = imap(i,1);
				u = imap(i,2);
				g = imap(k,1);
				h = imap(k,2);
				if v == g
					if u == h || u + 1 == h || u - 1 == h
						continue;
					end
				end

				dist = sqrt((y(t,x1_I)-y(t,x2_I))^2 + (y(t,y1_I)-y(t,y2_I))^2);
				LJval = LJatom_Energy(dist,epsi,sigma);
				Evdw = Evdw + LJval;
			end
		end
		%}

		Ebt(t) = Eb;
		Est(t) = Es;
		Evst(t) = Evs;
		Evdwt(t) = Evdw;
		E(t) = Eb + Es + Evs + Evdw;
	end
	
end

