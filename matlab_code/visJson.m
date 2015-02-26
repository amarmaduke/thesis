function visJson(json_file, sindex, varargin)
%visJson visualize simulation from JSON file
% visJson(file, sindex)
% visJson(file, sindex, [delay,figNumber,frameCount,radius])
%
% sindex is an assumed layer of complexity present in the
% json file that allows for several simulations to coexist
% in one file.
%
% delay is the wait time between updating to a new frame.
% if delay is -1 then a user-ended pause will occur.
%
% figNumber is the figure number attributed to the plot.
%
% frameCount is the number of frames to be iterated over.
% It is expected to be less than or equal to the maximum
% number of potential frames in the simulation data.
%
% radius will plot a circle of the desired radius around
% each point (for visualizing vdW cut-off).

p = inputParser;
defaultFigureNumber = 1;
defaultDelay = .1;
defaultRadius = 0;
defaultFrameCount = 100;

addRequired(p,'json_file');
addRequired(p,'sindex');
addOptional(p,'delay',defaultDelay);
addOptional(p,'figureNumber',defaultFigureNumber);
addOptional(p,'frameCount',defaultFrameCount);
addOptional(p,'radius',defaultRadius);

parse(p,json_file,sindex,varargin{:});
fign = p.Results.figureNumber;
frameCount = p.Results.frameCount;
pc = p.Results.delay;
radius = p.Results.radius;
json_file = p.Results.json_file;

state = struct();
json_struct = loadjson(json_file);
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

n = state.n;
m = state.m;
len = state.len;
N = state.n*state.m;
delta = state.delta
mvSubAtoms = 0:state.sub_h:(state.sub_count-1);

frameSkip = max(1,floor(length(y(:,1))/frameCount));

figure('units','normalized','outerposition',[0 0 1 1])
ma = N;
mn = -ma;

sa = [[-220 -150 0 10], [-150 -80 0 10], [-80 -10 0 10], [-10 60 0 10], [60 130 0 10]];
sa = [[-50, 50, 0 , 100]]

for	i = [1:frameSkip:length(y(:,1)), length(y(:,1))]
	sk = 1;
	for subk = 1:4:length(sa)
		subplot(length(sa)/4,1,sk);
		sk = sk + 1;
		hold on
		axis equal
		line([mn ma],[0 0],'Color','k','LineStyle','-');
		%axis equal
		%axis([mn ma -state.len (max(state.n)+1)*state.len]);
		%ma = 10; mn = -10;
		axis([sa(subk), sa(subk+1), sa(subk+2), sa(subk+3)]);

		xI = 2*N+1;
		yI = 2*N+2;
		s_pos(1) = y(i,xI);
		s_pos(2) = y(i,yI);
		plot(s_pos(1),s_pos(2),'o');
		for k = 1:length(mvSubAtoms)
			c_pos(1) = s_pos(1) + mvSubAtoms(k);
			c_pos(2) = s_pos(2);
			plot(c_pos(1),c_pos(2),'o');
			[c_x, c_y] = circ(c_pos(1),c_pos(2),radius);
			plot(c_x,c_y,':');
			%mn = min(mn,min(c_pos(1),c_pos(2)));
			%ma = max(ma,max(c_pos(1),c_pos(2)));
		end
		subMn = min(mvSubAtoms) + s_pos(1);
		subMa = max(mvSubAtoms) + s_pos(1);
		line([subMn subMa],[s_pos(2) s_pos(2)],'Color','k','LineStyle','-');

		c_pos(1) = state.osub;
		c_pos(2) = 0;
		for k = 1:state.osub_count
			plot(c_pos(1),c_pos(2),'o');
			[c_x, c_y] = circ(c_pos(1),c_pos(2),radius);
			plot(c_x,c_y,':');
			c_pos(1) = c_pos(1) + state.osub_h;
			c_pos(2) = c_pos(2);
		end

		for j = 1:m
			colorList = hsv(n);
			xI = 1 + (j-1)*n;
			yI = N + xI;
			c_pos(1) = delta(j);
			c_pos(2) = 0;
			plot(c_pos(1),c_pos(2),'o');
			for k = 1:n
				p_pos = c_pos;
				c_pos(1) = y(i,xI);
				c_pos(2) = y(i,yI);
				l_x = [ p_pos(1) c_pos(1) ]';
				l_y = [ p_pos(2) c_pos(2) ]';
				plot(c_pos(1),c_pos(2),'o','Color',colorList(k,:));
				[c_x, c_y] = circ(c_pos(1),c_pos(2),radius);
				plot(c_x,c_y,':','Color',colorList(k,:));
				mn = min(mn,min(c_pos(1),c_pos(2))+1);
				ma = max(ma,max(c_pos(1),c_pos(2))+1);
				plot(l_x,l_y,'Color','k');
				xI = xI + 1;
				yI = yI + 1;
			end
		end
	end
	
			%refreshdata
		%{
		if i == 1
			pause
		else
			if pc ~= -1
				pause(pc);
			else
				pause;
			end
		end
		%}
			
		drawnow
			
		frame = getframe(1);
		im = frame2im(frame);
		[imind, cm] = rgb2ind(im, 256);

		if i == 1
			imwrite(imind, cm, 'test.gif', 'gif', 'Loopcount', inf);
		else
			imwrite(imind, cm, 'test.gif', 'WriteMode', 'append');
		end

		%axis equal
		fprintf('Percent: %.2f Time: %.2f\n',((i/length(y(:,1))) * 100),(t(i)-1)*.1);

		if i ~= length(y(:,1))
			clf(fign);
		end
end

end