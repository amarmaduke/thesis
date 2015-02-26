function [ x, y ] = getPos( in, state, imap )
	
	delta = state.delta;
	N = state.N;
	cs = zeros(N,1);
	sn = zeros(N,1);
	for k = 1:length(cs)
		i = imap(k,2);
		if i == 1
			sn(k) = sin(in(k));
			cs(k) = cos(in(k));
		else
			sn(k) = sn(k-1) + sin(in(k));
			cs(k) = cs(k-1) + cos(in(k));
		end
	end
	
	x = zeros(N,1);
	y = zeros(N,1);

	for k = 1:N
		j = imap(k,1);
		x(k) = delta(j) + sn(k);
		y(k) = cs(k);
	end
	
end

