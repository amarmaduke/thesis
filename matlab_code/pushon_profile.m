linear_step = 5;
output_i = 1;
theta_counter = 1;
inner_counter = 0;
lower_bound = 0;
upper_bound = 100;
have_lower = false;
have_upper = false;
% 4:4:176
for theta = 4:4:176
	have_lower = false;
	have_upper = false;

	% 5:5:50
	for magnitude = 5:5:30
		done = false;
		load('sims/M2_seed.mat')
		inner_counter = 0;
		while ~done
			l = magnitude*sin(pi*theta/180);
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
			if equil < 1e-8 && adhesion < 1e-8 % Adhered
				[mvF, mvC, stF, stC] = compute_adhesion(y(end,:),state);
				done = true;
				output{output_i,1} = {theta, 1, l, m};
				output{output_i,2} = {theta, mvF, mvC, stF, stC};
				%output{output_i,1} = output_j;
				fprintf('Done\n');
			elseif abs(adhesion-d) < 1e-8  && inner_counter > 10 % Pulled Off
				[mvF, mvC, stF, stC] = compute_adhesion(y(end,:),state);
				done = true;
				output{output_i,1} = {theta, 0, l, m};
				output{output_i,2} = {theta, mvF, mvC, stF, stC};
				%output{output_i,1} = output_j;
				fprintf('Pulled Off\n');
			end
			inner_counter = inner_counter + 1;
		end
		
		temp = output{output_i,1};
		output_i = output_i + 1;
	end
	
	theta_counter = theta_counter + 1;
end

save('sims/M2_test','output');