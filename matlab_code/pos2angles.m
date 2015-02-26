function [out, lens] = pos2angles( u, state )

	delta = state.delta;
	m = state.m;
	n = state.n;
	N = state.N;
	
	out = zeros(1,N);
	lens = zeros(1,N);

	z = 1;
	for j = 1:m
		for i = 1:n(j)
			xI = z;
			yI = z + N;
			if i == 1
				dx = u(xI) - delta(j);
				dy = u(yI);
			else
				dx = u(xI) - u(xI - 1);
				dy = u(yI) - u(yI - 1);
			end
			theta = atan(dx / dy);
			if dx < 0 && dy < 0 || dy < 0
				theta = theta + pi;
			end
			l = sqrt(dx^2 + dy^2);
			out(z) = theta;
			lens(z) = l;
			
			z = z + 1;
		end
	end
end

