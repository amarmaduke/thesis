function F = genNghd( y, state, imap, radius, msIndex )
	
	r2 = radius*radius;
	N = state.N;
	mvSubAtoms = state.mvSubAtoms;
	stSubAtoms = state.stSubAtoms;
	stSubX = state.stSubX;
	n = state.n;
	
	k = 1;
	for i = 1:N
		for j = i+1:N
			v = imap(i,1);
			u = imap(i,2);
			g = imap(j,1);
			h = imap(j,2);
			% Rule out the atom itself and it's two immediate neighbors
			if v == g
				if u == h || u + 1 == h || u - 1 == h
					continue;
				end
			end
			
			if n(v) == 1 || n(g) == 1
				continue;
			end
			
			xps = y(i) - y(j);
			yps = y(i+N) - y(j+N);
			d = xps.^2 + yps.^2;
			
			if d <= r2
				C(k,1) = {i};
				C(k,2) = {j};
				k = k + 1;
			end
		end
	end
	
	if msIndex > 0
		mvSub_xI = 1 + sum(n(1:msIndex-1));
		mvSub_x = y(mvSub_xI);
		mvSub_y = y(mvSub_xI + N);
		for i = 1:N
			for j = 1:length(mvSubAtoms)
				if mvSub_xI == i
					continue;
				end
				
				xps = y(i) - (mvSub_x + mvSubAtoms(j));
				yps = y(i+N) - mvSub_y;
				
				d = xps.^2 + yps.^2;

				if d <= r2
					C(k,1) = {i};
					C(k,2) = {-j};
					k = k + 1;
				end
			end
		end
	end
	
	if ~isempty(stSubAtoms)
		for i = 1:N
			for j = 1:length(stSubAtoms)
				if mvSub_xI == i
					continue;
				end
				
				xps = y(i) - (stSubX + stSubAtoms(j));
				yps = y(i+N) - 0;
				
				d = xps.^2 + yps.^2;
				
				if d <= r2
					C(k,1) = {i};
					C(k,2) = {j*1i};
					k = k + 1;
				end
			end
		end
	end
	
	if exist('C','var')
		F = cell2mat(C);
	else
		F = [];
	end
end

