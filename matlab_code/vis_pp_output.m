M_adhered = 0; M = 0;
m_pulloff = 0; m = 1000000;
theta_v = 0;
reflect = false;
color = 'k';
fign = 6;

[K, ~] = size(output);

theta_p = -1;
theta_c = 0;
counter = 1;
for k = 1:K
	i = output{k,1};
	magnitude = sqrt(i{3}.^2 + i{4}.^2);
	theta_c = i{1};
	if k == 1
		theta_p = i{1};
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
	if i{2} == 1
		M = max(magnitude,M);
	else
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

%% Plot All
[K, ~] = size(output);
colorList = hsv(66);
figure(fign)
hold on
colormap hsv
for k = 1:K
	i = output{k,1};
	j = output{k,2};
	theta = i{1};
	mag = sqrt(i{3}^2 + i{4}^2);
	if i{2} == 0
		plot(theta,mag,'+','Color',colorList(j{3}+1,:))
	else
		plot(theta,mag,'x','Color',colorList(j{3}+1,:))
	end
end

ylabel('$\sqrt{\lambda^2 + \mu^2}$','interpreter','latex');
xlabel('$\theta$','interpreter','latex');

set(gca,'XDir','Reverse')

figure(fign+1)
hold on
for i = 1:33
	plot(i,i,'o','Color',colorList(i,:))
end

%% JSON Load
n_count = 96;
pstr = 'pullp_1_eb0';
str = ['Link to sims/n', num2str(n_count), '/runs/', pstr, '.json'];
%str = ['burnt.json'];
pp_struct = loadjson(str);
output = pp_struct.grid;

[sorted, sorted_order] = sort(output(:,1));
output = output(sorted_order,:);

fign = 4;

%% JSON

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

%% JSON plot All

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

%% JSON Plot lambda/mu grid
pp_struct = loadjson('sunout2.json');
output = pp_struct.grid;

[K, ~] = size(output);
colorList = hsv(66);
color = 'k';
figure(2)
hold on
colormap hsv
for k = 1:K
	i = output(k,:);
	if i(3) == 0
		plot(i(2),i(1),'+','Color',color)
	else
		plot(i(2),i(1),'x','Color',color)
	end
end

ylabel('$\lambda$','interpreter','latex');
xlabel('$\mu$','interpreter','latex');

set(gca,'YDir','Reverse')

%% JSON plot grid
pp_struct = loadjson('sunout2.json');
output = pp_struct.grid;

[K, ~] = size(output);
colorList = hsv(66);
figure(2)
hold on
colormap hsv
for k = 1:K
	i = output(k,:);
	theta = 180*atan2(i(2),-i(1))/pi + 90;
	mag = sqrt(i(1)^2 + i(2)^2);
	if i(3) == 0
		plot(theta,mag,'+','Color','k')
	else
		plot(theta,mag,'x','Color','k')
	end
end

ylabel('$\sqrt{\lambda^2 + \mu^2}$','interpreter','latex');
xlabel('$\theta$','interpreter','latex');

set(gca,'XDir','Reverse')

%% JSON plot Lambda/mu profile
pp_struct = loadjson('sunout2.json');
output = pp_struct.grid;

[K, ~] = size(output);
colorList = hsv(66);
color = 'k';
figure(4)
hold on
colormap hsv
for k = 1:K
	i = output(k,:);
	if i(4) == 0
		plot(i(3),i(2),'+','Color',color)
	else
		plot(i(3),i(2),'x','Color',color)
	end
end

ylabel('$\lambda$','interpreter','latex');
xlabel('$\mu$','interpreter','latex');

set(gca,'YDir','Reverse')

