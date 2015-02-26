function out = angles2pos( angles, lens, state )

	m = state.m;
	n = state.n;
	N = state.N;
	delta = state.delta;

	out = zeros(1,2*N);

	z = 1;
	for j = 1:m
		cs = 0;
		sn = delta(j);
		for i = 1:n(j)
			xI = z;
			yI = z + N;
			cs = cs + lens(z)*cos(angles(z));
			sn = sn + lens(z)*sin(angles(z));
			out(xI) = sn;
			out(yI) = cs;
			z = z + 1;
		end
	end

end

