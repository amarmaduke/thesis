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

%{
m = 3 + 1;
n = ones(1,m)*(m-1);
n(1) = 1;
beta = 1;
len = 1;
lambda = 0; mu = 0;
epsi = 1;
sigma = 1; 
gamma = 100;
delta = -1:(m-2);
density = 1;
mvSubAtoms = 0:30;
%}
r = 63.2; theta = 143.4;
state.lambda = -r*sin(theta*pi/180); state.mu = r*cos(theta*pi/180);
%state.lambda = -6.727-.1; state.mu = -12.87-.1;
mvSubIndex = 1;


checkInterval = 1000;
nhbdUpdate = 10;
nhbdRadius = 100000;

%%{
flag_IgnoreSetup = 0;
flag_GenerateInit = 0;
flag_ContinueFromPrevious = 1 & ~flag_GenerateInit;
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

N = sum(state.n);

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
	setup(1,1) = -5;
	setup(1,2) = 12; %m;
	%}
	init = [ setup(:,1); setup(:,2) ]';
end

if flag_ContinueFromPrevious && ~flag_IgnoreSetup
	[t_temp, s_temp] = size(y);
	init = y(t_temp,:);
end

if flag_RedoFromSave || flag_VisFromSave && ~flag_IgnoreSetup
	init = state.init;
end

%{
state = struct('m',m,'n',n,'N',N,'beta',beta,'len',len, ...
					'lambda',lambda,'mu',mu,'epsi',epsi,'sigma',sigma, ...
					'gamma',gamma,'delta',delta,'init',init, ...
					'mvSubAtoms',mvSubAtoms,'density',density);
%}
				
imap = zeros(state.N,2);
for i = 1:state.N
	[ imap(i,1), imap(i,2) ] = getAtom(i,state.n);
end

options = odeset('RelTol',1e-11,'AbsTol',1e-11);
%options = []

%}

%simpinit = [ init(1) init(2) init(5) init(6) ];

%% Solver and Time Domain
if ~flag_VisFromSave

state.tend = 100;
tspan = [ 0, state.tend-3, state.tend-2, state.tend-1, state.tend];
%tspan = 0:1:state.tend;
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
%vis(t,y,state,.1,1,10)
%vis(t2,y2,state,e2,.1)


%y = init;
%h = .0001;
%for i = 0:h:1
%	f = force(0,y,state,imap,checkInterval,nhbdUpdate,nhbdRadius,mvSubIndex);
%	y = y + h*f';
%end