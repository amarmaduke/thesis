%% Compute Using force.m
m = 2;
n = [1 1];
beta = 1;
len = 1;
lambda = 0; mu = 0;
epsi = 1;
sigma = 1; 
gamma = 100;
delta = [0 1];
density = 1;
mvSubAtoms = 0:2;

N = sum(n);

for j = 1:m
	for i = 1:n(j)
		setup(i + sum(n(1:j-1)),1) = delta(j);
		%if mod(i,2) == 0
		%	setup(i + sum(n(1:j-1)),1) = delta(j) - .25*rand(1,1);
		%end
		setup(i + sum(n(1:j-1)),2) = i;
	end
end
setup(1,1) = 0;
setup(1,2) = 1.9;
%}
init = [ setup(:,1); setup(:,2) ]';

state = struct('m',m,'n',n,'N',N,'beta',beta,'len',len, ...
					'lambda',lambda,'mu',mu,'epsi',epsi,'sigma',sigma, ...
					'gamma',gamma,'delta',delta,'init',init, ...
					'mvSubAtoms',mvSubAtoms,'density',density);

imap = zeros(state.N,2);
for i = 1:state.N
	[ imap(i,1), imap(i,2) ] = getAtom(i,state.n);
end

icheck = 1;
nUpdate = 1;
nRadius = 8;
msIndex = 1;
positions = init(1+N):.001:init(1+N)+3;
force_eval = zeros(1,length(positions));

for i = 1:length(positions)
	init(1+N) = positions(i);
	force_vector = forceEvaluation(0:1:length(init),init,state,imap,icheck,nUpdate,nRadius,msIndex);	
	force_eval(i) = sum(force_vector(4:5));
	%[ i / length(positions) ] 
end

hold on
figure(1)
plot(positions-1,force_eval,'o');

%% Compute Other Force
d = .9:.01:.9+3;
gamma = 1;
epsi = 1;
sigma = 1;

d1_1 = sqrt(.25*gamma^2 + d.^2);
F_1 = (24*epsi.*d./sigma./d1_1).*((sigma./d1_1).^13 - (sigma./d1_1).^7);

d1_2 = sqrt(gamma^2 + d.^2);
F_2 = (24*epsi.*d./sigma./d1_2).*((sigma./d1_2).^13 - (sigma./d1_2).^7) ...
		+ (12*epsi/sigma).*((sigma./d).^13 - (sigma./d).^7);

figure(2)
plot(d,F_1);

figure(1)
plot(d,F_2);

%% Timing Data

t = [3, 5, 7, 10, 15, 20, 25];
m_45 = [51.246, 178.139, 232.096, 534.418, 2231.347, 2832.672, nan];
m_15s = [3.492, 17.406, 28.022, 76.134, 628.769, 1682.726, 4917.949];
gpu = [34.864, ];

hold on
plot(t,m_45,'x',t,m_15s,'x');
plot(t,m_45,t,m_15s);

%% Stuff

output_i = 1;
for l = -25:0
	for m = -20:15
		output_j = 2;
		done = false;
		load('2014-2-11,10:40.mat')
		while ~done
			state.lambda = l; state.mu = m;
			main;

			[m, l]
			[output_i, output_j]
			output{output_i,output_j} = {t, y, state};
			output_j = output_j + 1;

			equil = norm(y(end,[2:state.N,state.N+2:end]) - ...
						 y(end-1,[2:state.N,state.N+2:end]))
			if equil < 1e-11 % Equillibrium
				adhesion = norm(y(end,[1,state.N+1])-y(end-1,[1,state.N+1]))
				d = norm([l,m])
				if abs(adhesion) < 1e-11 % Adhered
					done = true
					output{output_i,output_j} = {1, l, m};
					output{output_i,1} = output_j;
				elseif abs(adhesion - d) < 1e-11 % Free
					done = true
					output{output_i,output_j} = {0, l, m};
					output{output_i,1} = output_j;
				end
			end
		end
		output_i = output_i + 1;
	end
end

s = size(output);
figure(2)
hold on
for k = 1:s(1)
	last_j = output{k,1};
	C = output{k,last_j};
	x = C{3}; y = C{2};
	if C{1} == 0;
		plot(x,y,'x','Color','r');
	elseif C{1} == 1
		plot(x,y,'x','Color','k');
	end
end

axis equal
ylabel('$\lambda$.','interpreter','latex');
xlabel('$\mu$.','interpreter','latex');
axis equal
set(gca,'YDir','Reverse')
axis equal

%% Stuff

linear_step = 5;
output_i = 1;
theta_counter = 1;
magnitude = 75;
lower_bound = 0;
upper_bound = 100;
have_lower = false;
have_upper = false;
for theta = 0:2:20
	found = false;
	have_lower = false;
	have_upper = false;
	if theta ~= 0
		linear_step = 1;
	end
	while ~found
		output_j = 1;
		done = false;
		load('m1_seed.mat')
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
			if abs(adhesion - d) < 1e-11 % Pulled
				done = true;
				output{output_i,output_j} = {0, l, m};
				%output{output_i,1} = output_j;
				fprintf('Pulled Off\n');
			elseif equil < 1e-11 && adhesion < 1e-11 % Adhered
				done = true;
				output{output_i,output_j} = {1, l, m};
				%output{output_i,1} = output_j;
				fprintf('Adhered\n');
			end
		end
		
		temp = output{output_i,output_j};
		
		if have_upper && have_lower ...
			&& abs(upper_bound - lower_bound) ...
				/min(upper_bound,lower_bound) < .01
			found = true;
		end
		
		if ~have_lower || ~have_upper
			if temp{1} == 1
				lower_bound = magnitude;
				have_lower = true;
				if ~have_upper
					magnitude = magnitude + linear_step;
				end
			elseif temp{1} == 0
				upper_bound = magnitude;
				have_upper = true;
				if ~have_lower
					magnitude = magnitude - linear_step;
				end
			end
		end
		
		if have_lower && have_upper
			if temp{1} == 1
				lower_bound = magnitude;
			elseif temp{1} == 0
				upper_bound = magnitude;
			end
			magnitude = lower_bound + (upper_bound - lower_bound) / 2;
		end
		
		output_i = output_i + 1;
	end
	
	grid{theta_counter,1} = temp;
	%grid{theta_counter,2} = {output_i-1,output_j};
	theta_counter = theta_counter + 1;
end

s = size(output);
figure(10)
hold on
for k = 1:s(1)
	C = output{k};
	x = C{3}; y = C{2};
	if C{1} == 0;
		plot(x,y,'x','Color','r');
	elseif C{1} == 1
		plot(x,y,'x','Color','k');
	end
end

axis equal
ylabel('$\lambda$.','interpreter','latex');
xlabel('$\mu$.','interpreter','latex');
axis equal
set(gca,'YDir','Reverse')
axis equal

figure(11)
s = size(grid);
x = []; y = [];
theta = []; magnitude = [];
for k = 1:s(1)
	C = grid{k};
	x_0 = C{3}; y_0 = C{2};
	magnitude = [ magnitude, norm([x_0, y_0]) ];
	theta = [ theta, 180*atan2(-y_0,x_0)/pi ];
	x = [ x, x_0 ];
	y = [ y, y_0 ];
end

plot(x,y);

axis equal
ylabel('$\lambda$.','interpreter','latex');
xlabel('$\mu$.','interpreter','latex');
axis equal
set(gca,'YDir','Reverse')
axis equal

figure(12)
plot(theta,magnitude);
figure(2)
plot(theta,x);
figure(3)
plot(theta,y);
set(gca,'YDir','Reverse')
axis equal


%%
s = size(output);
figure(10)
hold on
for k = 1:s(1)
	last_j = output{k,1};
	C = output{k,last_j};
	x = C{3}; y = C{2};
	if C{1} == 0;
		plot(x,y,'x','Color','b');
	elseif C{1} == 1
		plot(x,y,'x','Color','g');
	end
end

axis equal
ylabel('$\lambda$.','interpreter','latex');
xlabel('$\mu$.','interpreter','latex');
axis equal
set(gca,'YDir','Reverse')
axis equal

figure(11)
s = size(grid);
x = []; y = [];
theta = []; magnitude = [];
for k = 1:s(1)
	C = grid{k,1};
	x_0 = C{3}; y_0 = C{2};
	magnitude = [ magnitude, norm([x_0, y_0]) ];
	theta = [ theta, 180*atan2(-y_0,x_0)/pi ];
	x = [ x, x_0 ];
	y = [ y, y_0 ];
end

plot(x,y);

axis equal
ylabel('$\lambda$.','interpreter','latex');
xlabel('$\mu$.','interpreter','latex');
axis equal
set(gca,'YDir','Reverse')
axis equal

figure(12)
plot(theta,magnitude);
figure(2)
plot(theta,x);
figure(3)
plot(theta,y);
set(gca,'YDir','Reverse')
axis equal

%% Parse JSON
json_file = 'temp1.json';
state = struct();
json_struct = loadjson(json_file);
fields = fieldnames(json_struct);
json_count = 1;
t = 0;
y = 0;
for i = 1:numel(fields)
	if strfind(fields{i},'tq')
		t(json_count) = sscanf(fields{i},'tq%d');
		y(json_count,:) = json_struct.(fields{i});
		json_count = json_count + 1;
	else
		state.(fields{i}) = json_struct.(fields{i});
	end
end

[t, t_order] = sort(t);
y = y(t_order,:);

%% Del

for k = 235:237
	output(235,:) = [];
end

%% Sort

t_array = cellfun(@(C) C{1}(1), output);
t_array = t_array(:,1);
[sorted, sorted_order] = sort(t_array);

output = output(sorted_order,:);

%%
r1 = 1;
rc = 3;
epsi = 100;
sigma = 1;

LJ = @(x) epsi*((sigma./x).^12 - 2*(sigma./x).^6);
LJp = @(x) (-12*epsi/sigma)*((sigma./x).^13 - (sigma./x).^7);
LJpp = @(x) (12*epsi/sigma/sigma)*(13*(sigma./x).^14 - 7*(sigma./x).^8);

U = @(x) LJ(x) - LJ(rc) - LJp(rc).*(x - rc);
Zero = @(x) 0;

F = @(x) -(LJp(x) - LJp(rc));



A = -LJ(rc);
B = LJp(rc);

t_alpha1 = (rc*rc - 2*r1*rc + r1*r1);
t_alpha2 = 1/(2*rc - 2*r1);
t_alpha3 = (3*rc*rc - 3*r1*r1);
t_alpha = -t_alpha1*t_alpha2*t_alpha3;

t_beta1 = (rc*rc - 2*rc*r1 +r1*r1);
t_beta2 = 1/(2*rc - 2*r1);
t_beta = t_beta1*t_beta2*B;

t_a = 1/(rc*rc*rc - 3*r1*r1*rc + 2*r1*r1*r1 - t_alpha);
a = t_a*(A - t_beta);

t_b1 = 1/(2*rc - 2*r1);
t_b2 = -(3*rc*rc - 3*r1*r1);
b = t_b1*(B - t_b2*a);

c = -(2*b*r1 + 3*a*r1*r1);

d = -(c*r1 + b*r1*r1 + a*r1*r1*r1);

P = @(x) a*x.^3 + b*x.^2 + c*x + d;

S = @(x) LJ(x) + P(x);

x = sigma-.2:.001:rc;
x2 = rc:.001:5;
y = P(x);
lj = LJ(x);
s = S(x);

y_f = LJ(x);
y_s = -LJp(x);
y_z = Zero(x2);

c = 'k';

figure(112);
hold on
plot(x,y_f,'Color',c);
line([rc, 5],[0, 0],'Color',c);
line([sigma-.2, 5],[-epsi, -epsi],'Color','k','LineStyle',':');
line([sigma, sigma],[-200, 1000],'Color','k','LineStyle',':');
axis([.8, 4, -110, 5]);

%line([subMn subMa],[s_pos(2) s_pos(2)],'Color','k','LineStyle','--');

figure(111);
hold on
plot(x,y_s,'Color',c);
line([rc, 5],[0, 0],'Color',c);
line([sigma, sigma],[-400, 1000],'Color','k','LineStyle',':');
axis([.8, 4, -300, 5]);

%figure(100)
%plot(x,y_f,'b',x,y_s,'r');
%plot(x,y,'b',x,lj,'k',x,s,'r');

%% Bioinfo Schtuff

alpha1 = 4;
epsi = 5;
beta = 6;
gamma = 40;
c = 1;
phi = 10000;

b1 = -2*c/epsi + (alpha1+2)*gamma/(epsi^(alpha1+1));
a1 = -gamma*(alpha1+1)/(epsi^(alpha1+2)) + c/(epsi^2);


P1 = @(x) a1*x.^(alpha1+2) + b1*x.^(alpha1+1) + c*x.^(alpha1);
E = @(x) gamma*exp(-phi*(x-epsi).^2);



x1 = -2*epsi:.001:epsi;
x2 = epsi:.001:2*epsi;

plot(x1,P1(x1),x2,E(x2));

%% Scripts
for j = 1:10
	for i = 1:96
		display([num2str(j*10), ','])
	end
end

for j = 1:10
	for i = 1:96
		display([num2str(i), ','])
	end
end

%% Mo money mo powa

sindex = 0;
state = struct();
json_struct = loadjson('burnt3.json');
fields = fieldnames(json_struct);
json_count = 1;
t2 = 0;
td = 0;
for i = 1:numel(fields)
	if strfind(fields{i},['sindex', num2str(sindex), 'tq'])
		t2(json_count) = sscanf(fields{i},['sindex', num2str(sindex), 'tq%d']);
		y2(json_count,:) = json_struct.(fields{i});
		json_count = json_count + 1;
	%elseif strfind(fields{i},['sindex', num2str(sindex), 'td'])
	%	td(json_count) = sscanf(fields{i},['sindex', num2str(sindex), 'td%d']);
	%	dy(json_count,:) = json_struct.(fields{i});
	%	json_count = json_count + 1;
	else
		state.(fields{i}) = json_struct.(fields{i});
	end
end

[t2, t_order] = sort(t2);
y2 = y2(t_order,:);

%% Plotuu~~
f1 = y(:,5);
f2 = y(:,5+96);
f3 = y(:,6);
f4 = y(:,6+96);
f5 = y(:,7);
f6 = y(:,7+96);
tz = t;
figure(4)
hold on
plot3(tz,f1,f2,'Color','b')
plot3(tz,f3,f4,'Color','r')
plot3(tz,f5,f6,'Color','g')
figure(5)
hold on
plot(tz,f1,'Color','b')
plot(tz,f3,'Color','r')
plot(tz,f5,'Color','g')
figure(6)
hold on
plot(tz,f2,'Color','b')
plot(tz,f4,'Color','r')
plot(tz,f6,'Color','g')

%% Fur loop
figure(2)
[t_size, p_size] = size(y);
[t_size2, p_size] = size(y2);
for i = 5000:t_size
	r = y(i,(1+96):(96+96));
	q = y(i,1:96);
	r2 = y2(i,(1+96):(96+96));
	q2 = y2(i,1:96);
	%q = 1:96;
	plot(q,r,'ro',q2,r2,'bx');
	axis([-100 100 0 3])
	pause(.01)
end
%%
[t_size, p_size] = size(y);
r = y(1:t_size,(1+96):(96+96));
q = y(1:t_size,1:96);
colormap(cool)
plot3(1:t_size,q,r,'x');
colormap(cool)
%%
for i = t_size+1:t_size2
	r2 = y2(i,(1+96):(96+96));
	q2 = y2(i,1:96);
	%q = 1:96;
	plot(q2,r2,'bx');
	axis([-100 100 0 3])
	pause(.01)
end
