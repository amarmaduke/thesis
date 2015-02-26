function out = forceDU(t,y,state,imap,icheck,nUpdate,nRadius,msIndex)
	persistent F;
	persistent DU;
	persistent fcount;
	persistent ltime;
	persistent ptime;

	n = state.n;
	m = state.m;
	beta = state.beta;
	len = state.len;
	lambda = state.lambda;
	mu = state.mu;
	epsi = state.epsi(1);
	epsi_top = state.epsi(3);
	epsi_bottom = state.epsi(2);
	sigma = state.sigma;
	gamma = state.gamma;
	delta = state.delta;
	density = state.density;
	N = state.N;
	mvSubAtoms = state.mvSubAtoms;
	stSubAtoms = state.stSubAtoms;
	stSubX = state.stSubX;
	
	if t == 0
		DU = y;
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
	
	% Check
	U = (DU(1:N) - y(1:N)).^2 + (DU(N+1:end) - y(N+1:end)).^2;
	u_ = max(sqrt(U));
	
	if u_ >= nRadius/2 || fcount == 1
		DU = y;
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
			
			b_3_t1 = 4*(lens(l_i+2)/lens(l_i+1))*product_f;
			b_3_t2 = 4*(lens(l_i+2)*lens(l_i+1));
			b_3_b_t = lens(l_i+2)*lens(l_i+1) + product_f;
			b_3_b = b_3_b_t*b_3_b_t;
			b_3x = (b_3_t1*(pvec(xI)-pvec(xI+1)) - ... 
					b_3_t2*(pvec(xI+1)-pvec(xI+2)))/b_3_b;
			b_3y = (b_3_t1*(pvec(yI)-pvec(yI+1)) - ... 
					b_3_t2*(pvec(yI+1)-pvec(yI+2)))/b_3_b;

			b_2_t1x = 4*(lens(l_i)/lens(l_i+1)*(pvec(xI)-pvec(xI+1)) + ...
				lens(l_i+1)/lens(l_i)*(pvec(xI)-pvec(xI-1)))*product_c;
			b_2_t2x = 4*lens(l_i+1)*lens(l_i)* ...
						(pvec(xI-1)-2*pvec(xI)+pvec(xI+1));
			b_2_t1y = 4*(lens(l_i)/lens(l_i+1)*(pvec(yI)-pvec(yI+1)) + ...
				lens(l_i+1)/lens(l_i)*(pvec(yI)-pvec(yI-1)))*product_c;
			b_2_t2y = 4*lens(l_i+1)*lens(l_i)* ...
						(pvec(yI-1)-2*pvec(yI)+pvec(yI+1));
			b_2_b_t = lens(l_i+1)*lens(l_i) + product_c;
			b_2_b = b_2_b_t*b_2_b_t;
			b_2x = (b_2_t1x - b_2_t2x)/b_2_b;
			b_2y = (b_2_t1y - b_2_t2y)/b_2_b;
			
			b_1_t1 = 4*(lens(l_i-1)/lens(l_i))*product_b;
			b_1_t2 = 4*(lens(l_i)*lens(l_i-1));
			b_1_b_t = lens(l_i)*lens(l_i-1) + product_b;
			b_1_b = b_1_b_t*b_1_b_t;
			b_1x = (b_1_t1*(pvec(xI)-pvec(xI-1)) - ... 
					b_1_t2*(pvec(xI-1)-pvec(xI-2)))/b_1_b;
			b_1y = (b_1_t1*(pvec(yI)-pvec(yI-1)) - ... 
					b_1_t2*(pvec(yI-1)-pvec(yI-2)))/b_1_b;
			
			if isnan(b_1x), b_1x = 0; end
			if isnan(b_2x), b_2x = 0; end
			if isnan(b_3x), b_3x = 0; end
			
			if isnan(b_1y), b_1y = 0; end
			if isnan(b_2y), b_2y = 0; end
			if isnan(b_3y), b_3y = 0; end
			
			xout = b_1x + b_2x + b_3x;
			yout = b_1y + b_2y + b_3y;
			
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
	
	%% V.W. atom-to-atom and atom-to-substrate
	%%{
	for h = 1:length(F)
		j = F(h,1);
		i = F(h,2);

		if imag(i) ~= 0
			i = imag(i);
			
			xps = stSubX + stSubAtoms(i) - y(j);
			yps = 0 - y(j+N);
			
			d = xps.^2 + yps.^2;
			dist = sqrt(d);
			temp_x = xps/dist;
			temp_y = yps/dist;
			LJval = LJatom(dist,epsi_bottom,sigma);
			
			out(j) = out(j) + LJval*temp_x;
			out(j+N) = out(j+N) + LJval*temp_y;
		elseif i < 0
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
			LJval = LJatom(dist,epsi_top,sigma);
				
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