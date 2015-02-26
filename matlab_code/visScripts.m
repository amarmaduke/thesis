%% JSON Load (Pulloff/Pushon Profile)
n_count = 96;
pstr = 'pullp3_3';
str = ['Link to sims/n', num2str(n_count), '/runs/', pstr, '.json'];
%str = ['temp.json'];
pp_struct = loadjson(str);
output = pp_struct.grid;

[sorted, sorted_order] = sort(output(:,1));
output = output(sorted_order,:);

fign = 4;

%% JSON Visualize Pulloff Boundary

M_adhered = 0; M = 0;
m_pulloff = 0; m = 1000000;
theta_v = 0;
color = 'k';
reflect = false;

[K, ~] = size(output);

theta_p = -1;
theta_c = 0;
counter = 1;
for k = 1:K
	i = output(k,:);
	magnitude = sqrt(i(2).^2 + i(3).^2);
	theta_c = i(1);
	if k == 1
		theta_p = i(1);
	end
	
	if theta_c ~= theta_p || k == K
		M_adhered(counter) = M;
		m_pulloff(counter) = m;
		if reflect && theta_p > 90
			theta_v(counter) = 180 - theta_p;
		else
			theta_v(counter) = theta_p;
		end
		counter = counter + 1;
		M = 0;
		m = 1000000;
		theta_p = theta_c;
	end
	
	if i(4) == 1
		M = max(magnitude,M);
	elseif i(4) == 0
		m = min(magnitude,m);
	end
end

figure(fign)
hold on
plot(theta_v,M_adhered,'--x','Color',color)
plot(theta_v,m_pulloff,color)
plot(theta_v,m_pulloff,'+','Color',color)

ylabel('$\sqrt{\lambda^2 + \mu^2}$','interpreter','latex');
xlabel('$\theta$','interpreter','latex');

set(gca,'XDir','Reverse')

%% JSON plot All (Pulloff or Pushon)

[K, ~] = size(output);
colorList = cool(n_count+1);
color = 'k';
figure(fign)
hold on
h = gcf;
dc_obj = datacursormode(h);
set(dc_obj,'UpdateFcn',@(obj, eobj) updateFcn(output, eobj))
colormap(colorList)
for k = 1:K
	i = output(k,:);
	theta = i(1);
	mag = sqrt(i(2)^2 + i(3)^2);
	if i(4) == 0
		plot(theta,mag,'+','Color',colorList(i(7)+1,:),'UserData',k);
	elseif i(4) == 1
		plot(theta,mag,'x','Color',colorList(i(7)+1,:),'UserData',k)
	elseif i(4) == 2
		plot(theta,mag,'o','Color',colorList(i(7)+1,:),'UserData',k)
	end
end

ylabel('$\sqrt{\lambda^2 + \mu^2}$','interpreter','latex');
xlabel('$\theta$','interpreter','latex');

set(gca,'XDir','Reverse')

%% JSON plot all (free-standing)
[K, ~] = size(output);
colorList = cool(n_count+1);
color = 'k';
figure(fign)
hold on
h = gcf;
dc_obj = datacursormode(h);
set(dc_obj,'UpdateFcn',@(obj, eobj) updateFcn(output, eobj))
colormap(colorList)
for k = 1:K
	i = output(k,:);
	b = i(1);
	eb = i(2);
	plot(eb,b,'.','Color',colorList(i(6)+1,:),'UserData',k);
end

%% Plotters we used before to visualize the waves
%% Load Json
sindex = 0;
state = struct();
json_struct = loadjson('0.json');
fields = fieldnames(json_struct);
json_count = 1;
t = 0;
td = 0;
for i = 1:numel(fields)
	if strfind(fields{i},['sindex', num2str(sindex), 'tq'])
		t(json_count) = sscanf(fields{i},['sindex', num2str(sindex), 'tq%d']);
		y(json_count,:) = json_struct.(fields{i});
		json_count = json_count + 1;
	%elseif strfind(fields{i},['sindex', num2str(sindex), 'td'])
	%	td(json_count) = sscanf(fields{i},['sindex', num2str(sindex), 'td%d']);
	%	dy(json_count,:) = json_struct.(fields{i});
	%	json_count = json_count + 1;
	else
		state.(fields{i}) = json_struct.(fields{i});
	end
end

[t, t_order] = sort(t);
y = y(t_order,:);
%% Individual points through time

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

%% Movie of the fiber through time
figure(2)
[t_size, p_size] = size(y);
%[t_size2, p_size] = size(y2);
for i = 1:t_size
	r = y(i,(1+96):(96+96));
	q = y(i,1:96);
%	r2 = y2(i,(1+96):(96+96));
%	q2 = y2(i,1:96);
	%q = 1:96;
	plot(q,r,'ro');
	%plot(q,r,'ro',q2,r2,'bx');
	axis([-100 100 0.25 2])
	pause(.01)
end

%{
for i = t_size+1:t_size2
	r2 = y2(i,(1+96):(96+96));
	q2 = y2(i,1:96);
	%q = 1:96;
	plot(q2,r2,'bx');
	axis([-100 100 0 3])
	pause(.01)
end
%}

%% 3D (x,y,time)

colorList = hsv(96+1);
[t_size, p_size] = size(y);
%r = y(1:t_size,(1+96):(96+96));
%q = y(1:t_size,1:96);
t = 1:t_size;
hold on
for i = 1:96
	x1 = y(:,i);
	y1 = y(:,i+96);
	plot3(t,x1,y1,'Marker','x','Color',colorList(i,:))
end