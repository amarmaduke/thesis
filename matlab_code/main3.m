%% TODO:

% Atomize the bottom substrate

%% Mapping
%	Given (j, i) ~ ( tube index, atom index )
%	and an array of positions, y, then
%	y(i + sum(n(1:j-1))) is the x position for the atom
%	at (j, i) and
%	y(i + sum(n(1:j-1)) + sum(n)) is the y position
%% Parameters
%	m := # of tubes
%	n(j) := # of atoms on j-th tube
%	beta := bending coefficient
%	len := desired length between atoms on a tube
%	gamma := strenght of spring between atoms
%	lambda := vertical load on the moving substrate
%	mu := horizontal load on the moving substrate
%	epsi := Strength of Van der Waals interaction
%	sigma := Equillibrium distance of Van der Waals.
%	init := initial condition
%	delta(j) := offset of the j-th tube
%	density := LJ substrate atom density
%	mvSubAtoms := position of atoms on the moving substrate

%%{
m = 1 + 1;
n = ones(1,m)*96;%(m-1);
n(1) = 1;
beta = 1;
len = 1;
%r = 4.613+.1; theta = 83.6;
%lambda = -r*sin(theta*pi/180); mu = r*cos(theta*pi/180);
%lambda = -6.727-.1; mu = -12.87-.1;
lambda = 10; mu = 0;
epsi_bottom = 1;
epsi_top = 1;
epsi = [1, epsi_bottom, epsi_top];
sigma = 1;
gamma = 100;
delta = [1, 0];
density = 0;
mvSubAtoms = 0:149;
stSubAtoms = 0:149;
stSubX = -75;
%}
mvSubIndex = 1;


checkInterval = 1000;
nhbdUpdate = 1000;
nhbdRadius = 8;

%%{
flag_IgnoreSetup = 0;
flag_GenerateInit = 1;
flag_ContinueFromPrevious = 0 & ~flag_GenerateInit;
flag_RedoFromSave = 0;
flag_VisFromSave = 0;
%}

%{
if length(delta) ~= length(n) || length(n) ~= m || m ~= length(delta)
	error('delta, n, or m incorrect');
end

if isempty(find(n == 1, 1)) && mvSubIndex < 0
	warninig('You seem to want a moving substrate but have not specified its index')
end
%}

N = sum(n);

if flag_GenerateInit && ~flag_IgnoreSetup
	clear('setup');

	%%{
	rng(1)
	for j = 1:m
		for i = 1:n(j)
			setup(i + sum(n(1:j-1)),1) = delta(j);
			%if mod(i,2) == 0
			%	setup(i + sum(n(1:j-1)),1) = delta(j) - .25*rand(1,1);
			%end
			setup(i + sum(n(1:j-1)),2) = i;
		end
	end
	setup(1,1) = -75;
	setup(1,2) = 100;
	%}
	init = [ setup(:,1); setup(:,2) ]';
end

%init = [sans.init(193), sans.init(1:96), sans.init(194), sans.init(97:192)];
%per = .0000001;
%init = [-5, 1:17, 15:-1:1, 33, 1:17, ones(1,15)*17];
%init_o = [1:17, 15:-1:1, 1:17, ones(1,15)*17, -5, 33];

if flag_ContinueFromPrevious && ~flag_IgnoreSetup
	[t_temp, s_temp] = size(y);
	init = y(t_temp,:);
end


if flag_RedoFromSave || flag_VisFromSave && ~flag_IgnoreSetup
	init = state.init;
end

%%{
state = struct('m',m,'n',n,'N',N,'beta',beta,'len',len, ...
					'lambda',lambda,'mu',mu,'epsi',epsi,'sigma',sigma, ...
					'gamma',gamma,'delta',delta,'init',init, ...
					'mvSubAtoms',mvSubAtoms,'density',density, ...
					'stSubAtoms',stSubAtoms,'stSubX',stSubX);
%}
				
imap = zeros(state.N,2);
for i = 1:state.N
	[ imap(i,1), imap(i,2) ] = getAtom(i,state.n);
end

options = odeset('RelTol',1e-12,'AbsTol',1e-12);
%options = []

%}

%{


s_ = n(2)*(m-1);
init = init_t([2*s_+1,1:s_,2*s_+2,s_+1:2*s_]);
%simpinit = [ init(1) init(2) init(5) init(6) ];
f_eval = force(0,init,state,imap,checkInterval,nhbdUpdate,nhbdRadius,mvSubIndex);
f = f_eval([2:s_+1,s_+3:2*s_+2,1,s_+2]);
%}
%{
s_ = n(2)*(m-1);
init_t = [-9.0622896e-21,-5.1939698e-20,-7.9668641e-19,-5.4560089e-17,6.2221462e-17,-5.5487796e-17,4.7139067e-17,1.5981336e-18,-4.9293089e-20,-5.5912607e-22,1.0346524,1.9488216,3.2521307,3.6534085,5.5011078,5.5020527,7.4760441,7.7162316,9.1649077,9.9649873,-100,100];
init = init_t([2*s_+1,1:s_,2*s_+2,s_+1:2*s_]);
f_eval = forceDU(0,init,state,imap,checkInterval,nhbdUpdate,nhbdRadius,mvSubIndex);
f = f_eval([2:s_+1,s_+3:2*s_+2,1,s_+2]);
%}

%{
s_ = n(2)*(m-1);
f_eval = forceDU(0,init,state,imap,checkInterval,nhbdUpdate,nhbdRadius,mvSubIndex);
f = f_eval([2:s_+1,s_+3:2*s_+2,1,s_+2]);
%}

%% Solver and Time Domain
%%{
if ~flag_VisFromSave

state.tend = 3;
%tspan = [ 0, state.tend-3, state.tend-2, state.tend-1, state.tend];
tspan = 0:1:state.tend;
startTime = cputime;
%[ t2, y2 ] = ode45(@simpForce,tspan,init,options,state);
[t, y] = ode15s(@forceDU,tspan,init,options,state,imap, ... 
					checkInterval,nhbdUpdate,nhbdRadius,mvSubIndex);
%[t2, y2] = ode45(@force,tspan,init,options,state,imap, ... 
%					checkInterval,nhbdUpdate,nhbdRadius,mvSubIndex);

endTime = cputime;
elapsedTime = endTime - startTime;
state.time = elapsedTime;


%for i = 1:state.tend
%	err = max(y(i,:)-y2(i,:));
%	if err > 1e-6
%		[i, err]
%	end
%end

end

%[e, b, s, sv, v] = energy(y,state,imap);
%[e2, b, s, sv, v] = energy(y2,state,imap);
%[e2, b2, s2, sv2] = simpEnergy(y,state);
%close all
%plot(t,e,'k',t,b,'r',t,s,'b',t,sv,':g',t,v,'g');
%plot(t,e);
vis(t,y,state,.1,2,10)
%vis(t2,y2,state,e2,.1)


%y = init;
%h = .0001;
%for i = 0:h:1
%	f = force(0,init,state,imap,checkInterval,nhbdUpdate,nhbdRadius,mvSubIndex);
%	y = y + h*f';
%end
%}