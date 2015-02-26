function out = force(t,y,state,imap,icheck,nUpdate,nRadius,msIndex)
	persistent F;
	persistent fcount;
	persistent ltime;
	persistent ptime;

	n = state.n;
	m = state.m;
	beta = state.beta;
	len = state.len;
	lambda = state.lambda;
	mu = state.mu;
	epsi = state.epsi;
	sigma = state.sigma;
	gamma = state.gamma;
	delta = state.delta;
	density = state.density;
	N = state.N;
	mvSubAtoms = state.mvSubAtoms;
	
	if t == 0
		ltime = cputime;
		fcount = 0;
		ptime = 0;
	end
	fcount = fcount + 1;
	
	%% Timing
	if mod(fcount,icheck) == 0
		ttime = cputime;
		elpsd = ttime - ltime;
		ptime = ptime + elpsd;
		pcnt = 100 * (t / state.tend);
		cmpltnTime = ( ( ptime * 100 / pcnt ) - ptime ) / 60;
		fprintf('%.2f %%, %.5f s, %.5f m; %.5f m \n', pcnt, elpsd, ...
			ptime/60,cmpltnTime);
		ltime = cputime;
	end
	
	%% Update Neighborhood
	%%{
	if mod(fcount,nUpdate) == 0 || fcount == 1
		F = genNghd(y,state,imap,nRadius,msIndex);
	end
	%}
	
	out = zeros(length(y),1);
	
	%% Pad Vector
	%%{
	pvec = zeros(2*(N+4),1);
	px_start = zeros(m,1);
	py_start = zeros(m,1);
	z = 1;
	for j = 1:m % Setup x padding
		px_start(j) = z+2;
		pvec(z) = delta(j);
		pvec(z+1) = delta(j);
		y_start = 1 + sum(n(1:j-1));
		y_end = y_start + n(j) - 1;
		pvec(z+2:z+n(j)+1) = y(y_start:y_end);
		pvec(z+n(j)+2) = NaN;
		pvec(z+n(j)+3) = NaN;
		z = z + n(j) + 4;
	end
	
	for j = 1:m % Setup y padding
		py_start(j) = z+2;
		pvec(z) = -len;
		pvec(z+1) = 0;
		y_start = 1 + sum(n(1:j-1)) + N;
		y_end = y_start + n(j) - 1;
		pvec(z+2:z+n(j)+1) = y(y_start:y_end);
		pvec(z+n(j)+2) = NaN;
		pvec(z+n(j)+3) = NaN;
		z = z + n(j) + 4;
	end
	%}
	
	%% Compute Lengths
	z = 1;
	lens = zeros(N+3*m,1);
	lens_start = zeros(m,1);
	for j = 1:m
		xI = -1 + px_start(j); % x index for vector (j, i)
		yI = -1 + py_start(j); % y index
			
		lens(z) = sqrt((pvec(xI) - pvec(xI-1))^2 + ...
						(pvec(yI) - pvec(yI-1))^2);
		z = z + 1;
		lens_start(j) = z;
		for i = 1:n(j)
			xI = (i-1) + px_start(j); % x index for vector (j, i)
			yI = (i-1) + py_start(j); % y index
			
			lens(z) = sqrt((pvec(xI) - pvec(xI-1))^2 + ...
						(pvec(yI) - pvec(yI-1))^2);
			
			z = z + 1;
		end
		lens(z) = NaN; lens(z+1) = NaN;
		z = z + 2;
	end
	
	%% V.W. /w Substrate + Bending + Load + Spring
	%%{
	z = 1;
	for j = 1:m
		for i = 1:n(j)
			xI = (i-1) + px_start(j); % x index for vector (j, i)
			yI = (i-1) + py_start(j); % y index
			xI_out = i + sum(n(1:j-1));
			yI_out = xI_out + N;
			l_i = (i-1) + lens_start(j);
				
			% Load
			%%{
			if n(j) == 1
				out(xI_out) = out(xI_out) + mu;
				out(yI_out) = out(yI_out) - lambda;
				z = z + 1;
				continue;
			end
			%}
			
			% Bending
			%%{
			xd_f = (pvec(xI+1) - pvec(xI))*(pvec(xI+2) - pvec(xI+1));
			yd_f = (pvec(yI+1) - pvec(yI))*(pvec(yI+2) - pvec(yI+1));
			xd_c = (pvec(xI) - pvec(xI-1))*(pvec(xI+1) - pvec(xI));
			yd_c = (pvec(yI) - pvec(yI-1))*(pvec(yI+1) - pvec(yI));
			xd_b = (pvec(xI-1) - pvec(xI-2))*(pvec(xI) - pvec(xI-1));
			yd_b = (pvec(yI-1) - pvec(yI-2))*(pvec(yI) - pvec(yI-1));
			
			product_f = xd_f + yd_f;
			product_c = xd_c + yd_c;
			product_b = xd_b + yd_b;
			
			forward_1 = lens(l_i+2) / (lens(l_i+1)*product_f);
			forward_2 = lens(l_i+1)*lens(l_i+2) / (product_f^2);
			center_1 = lens(l_i+1) / (lens(l_i)*product_c);
			center_2 = lens(l_i) / (lens(l_i+1)*product_c);
			center_3 = -lens(l_i)*lens(l_i+1) / (product_c^2);
			backward_1 = lens(l_i-1) / (lens(l_i)*product_b);
			backward_2 = lens(l_i-1)*lens(l_i) / (product_b^2);
			
			forward_1x = forward_1*(pvec(xI) - pvec(xI+1));
			forward_2x = forward_2*(pvec(xI+2) - pvec(xI+1));
			center_1x = center_1*(pvec(xI) - pvec(xI-1));
			center_2x = center_2*(pvec(xI) - pvec(xI+1));
			center_3x = center_3*(pvec(xI+1) - 2*pvec(xI) + pvec(xI-1));
			backward_1x = backward_1*(pvec(xI) - pvec(xI-1));
			backward_2x = backward_2*(pvec(xI-2) - pvec(xI-1));
			
			forward_1y = forward_1*(pvec(yI) - pvec(yI+1));
			forward_2y = forward_2*(pvec(yI+2) - pvec(yI+1));
			center_1y = center_1*(pvec(yI) - pvec(yI-1));
			center_2y = center_2*(pvec(yI) - pvec(yI+1));
			center_3y = center_3*(pvec(yI+1) - 2*pvec(yI) + pvec(yI-1));
			backward_1y = backward_1*(pvec(yI) - pvec(yI-1));
			backward_2y = backward_2*(pvec(yI-2) - pvec(yI-1));
			
			if isnan(forward_1x), forward_1x = 0; end
			if isnan(forward_2x), forward_2x = 0; end
			if isnan(center_1x), center_1x = 0; end
			if isnan(center_2x), center_2x = 0; end
			if isnan(center_3x), center_3x = 0; end
			if isnan(backward_1x), backward_1x = 0; end
			if isnan(backward_2x), backward_2x = 0; end
			
			if isnan(forward_1y), forward_1y = 0; end
			if isnan(forward_2y), forward_2y = 0; end
			if isnan(center_1y), center_1y = 0; end
			if isnan(center_2y), center_2y = 0; end
			if isnan(center_3y), center_3y = 0; end
			if isnan(backward_1y), backward_1y = 0; end
			if isnan(backward_2y), backward_2y = 0; end
			
			forward_x = forward_1x + forward_2x;
			center_x = center_1x + center_2x + center_3x;
			backward_x = backward_1x + backward_2x;
					
			forward_y = forward_1y + forward_2y;
			center_y = center_1y + center_2y + center_3y;
			backward_y = backward_1y + backward_2y;
			
			xout = forward_x + center_x + backward_x;
			yout = forward_y + center_y + backward_y;
			
			out(xI_out) = out(xI_out) - beta*xout;
			out(yI_out) = out(yI_out) - beta*yout;
			%}
			
			% Spring
			%%{
			forward = (lens(l_i+1) - len)/lens(l_i+1);
			backward = (lens(l_i) - len)/lens(l_i);
			
			forward_x = forward*2*(pvec(xI) - pvec(xI+1));
			backward_x = backward*2*(pvec(xI) - pvec(xI-1));
			if isnan(forward_x), forward_x = 0; end
			if isnan(backward_x), backward_x = 0; end
			
			forward_y = forward*2*(pvec(yI) - pvec(yI+1));
			backward_y = backward*2*(pvec(yI) - pvec(yI-1));
			if isnan(forward_y), forward_y = 0; end
			if isnan(backward_y), backward_y = 0; end
			
			out(xI_out) = out(xI_out) - gamma*(forward_x + backward_x);
			out(yI_out) = out(yI_out) - gamma*(forward_y + backward_y);
			%}
			
			% V.W. /w Substrate
			out(yI_out) = out(yI_out) - LJsubstrate(y(yI_out),epsi,sigma,density);
			
			z = z + 1;
		end
	end
	
	%% V.W. atom-to-atom
	%%{
	for h = 1:length(F)
		j = F(h,1);
		i = F(h,2);
		
		if i < 0
			i = abs(i);
			mvSub_xI = 1 + sum(n(1:msIndex-1));
			mvSub_yI = mvSub_xI + N;
			mvSub_x = y(mvSub_xI);
			mvSub_y = y(mvSub_yI);
			
			xps = mvSub_x + mvSubAtoms(i) - y(j);
			yps = mvSub_y - y(j+N);
			d = xps.^2 + yps.^2;
			dist = sqrt(d);
			temp_x = xps/dist;
			temp_y = yps/dist;
			LJval = LJatom(dist,epsi,sigma);
				
			out(mvSub_xI) = out(mvSub_xI) - LJval*temp_x;
			out(mvSub_yI) = out(mvSub_yI) - LJval*temp_y;
			out(j) = out(j) + LJval*temp_x;
			out(j+N) = out(j+N) + LJval*temp_y;
		else
			xps = y(i) - y(j);
			yps = y(i+N) - y(j+N);
			d = xps.^2 + yps.^2;
			dist = sqrt(d);
			temp_x = xps/dist;
			temp_y = yps/dist;
			LJval = LJatom(dist,epsi,sigma);
			
			out(i) = out(i) - LJval*temp_x;
			out(i+N) = out(i+N) - LJval*temp_y;
			out(j) = out(j) + LJval*temp_x;
			out(j+N) = out(j+N) + LJval*temp_y;
		end
	end
	%}
end