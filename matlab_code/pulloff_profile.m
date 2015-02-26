linear_step = 10;
output_i = 1;
theta_counter = 1;
magnitude = 65;
lower_bound = 0;
upper_bound = 100;
have_lower = false;
have_upper = false;
for theta = 140:4:160
	found = false;
	have_lower = false;
	have_upper = false;

	while ~found
		done = false;
		load('sims/m1_seed.mat')
		while ~done
			l = -magnitude*sin(pi*theta/180);
			m = magnitude*cos(pi*theta/180);
			state.lambda = l; state.mu = m;

			fprintf('------------------------------\n');
			fprintf('theta: %d, mag: %f\n',theta,magnitude)
			fprintf('m: %f, l: %f\n',m,l)
			fprintf('------------------------------\n');
			
			main;

			%output{output_i,output_j} = {t, y, state};
			%output_j = output_j + 1;

			equil = norm(y(end,[2:state.N,state.N+2:end]) - ...
						 y(end-1,[2:state.N,state.N+2:end]));
			adhesion = norm(y(end,[1,state.N+1])-y(end-1,[1,state.N+1]));
			d = norm([l,m]);
			if abs(adhesion - d) < 1e-8 % Pulled
				[mvF, mvC, stF, stC] = compute_adhesion(y(end,:),state);
				done = true;
				output{output_i,1} = {theta, 0, l, m};
				output{output_i,2} = {theta, mvF, mvC, stF, stC};
				%output{output_i,1} = output_j;
				fprintf('Pulled Off\n');
			elseif equil < 1e-11 && adhesion < 1e-8 % Adhered
				[mvF, mvC, stF, stC] = compute_adhesion(y(end,:),state);
				done = true;
				output{output_i,1} = {theta, 1, l, m};
				output{output_i,2} = {theta, mvF, mvC, stF, stC};
				%output{output_i,1} = output_j;
				fprintf('Adhered\n');
			end
		end
		
		temp = output{output_i,1};
		
		if have_upper && have_lower ...
			&& abs(upper_bound - lower_bound) ...
				/min(upper_bound,lower_bound) < .01
			found = true;
		end
		
		if ~have_lower || ~have_upper
			if temp{2} == 1
				lower_bound = magnitude;
				have_lower = true;
				if ~have_upper
					magnitude = magnitude + linear_step;
				end
			elseif temp{2} == 0
				upper_bound = magnitude;
				have_upper = true;
				if ~have_lower
					magnitude = magnitude - linear_step;
				end
			end
		end
		
		if have_lower && have_upper
			if temp{2} == 1
				lower_bound = magnitude;
			elseif temp{2} == 0
				upper_bound = magnitude;
			end
			magnitude = lower_bound + (upper_bound - lower_bound) / 2;
		end
		
		output_i = output_i + 1;
	end
	
	theta_counter = theta_counter + 1;
end

save('m1_test','output');