function vis(t, y, state, varargin)

p = inputParser;
defaultFigureNumber = 1;
defaultDelay = .1;
defaultRadius = 0;
defaultFrameCount = 100;

addRequired(p,'t');
addRequired(p,'y');
addRequired(p,'state');
addOptional(p,'delay',defaultDelay);
addOptional(p,'figureNumber',defaultFigureNumber);
addOptional(p,'frameCount',defaultFrameCount);
addOptional(p,'radius',defaultRadius);

parse(p,t,y,state,varargin{:});
fign = p.Results.figureNumber;
frameCount = p.Results.frameCount;
pc = p.Results.delay;
radius = p.Results.radius;
t = p.Results.t;
y = p.Results.y;
state = p.Results.state;

n = state.n;
len = state.len;
N = state.N;
delta = state.delta;
mvSubAtoms = state.mvSubAtoms;

frameSkip = max(1,floor(length(y(:,1))/frameCount));

figure(fign)
ma = N;
mn = -ma;
for i = 1:frameSkip:length(y(:,1))
	hold on
	axis equal
	line([mn ma],[0 0],'Color','k','LineStyle','-');
	%axis equal
	%axis([mn ma -state.len (max(state.n)+1)*state.len]);
	%ma = 10; mn = -10;
	axis([mn ma 0 ma]);
	
	for j = 1:length(state.n)
		colorList = hsv(n(j));
		if n(j) == 1
			xI = 1 + sum(n(1:j-1));
			yI = N + xI;
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
		else
			xI = 1 + sum(n(1:j-1));
			yI = N + xI;
			c_pos(1) = delta(j);
			c_pos(2) = 0;
			plot(c_pos(1),c_pos(2),'o');
			for k = 1:state.n(j)
				p_pos = c_pos;
				c_pos(1) = y(i,xI);
				c_pos(2) = y(i,yI);
				l_x = [ p_pos(1) c_pos(1) ]';
				l_y = [ p_pos(2) c_pos(2) ]';
				plot(c_pos(1),c_pos(2),'o','Color',colorList(k,:));
				[c_x, c_y] = circ(c_pos(1),c_pos(2),radius);
				plot(c_x,c_y,':','Color',colorList(k,:));
				mn = min(mn,min(c_pos(1),c_pos(2)));
				ma = max(ma,max(c_pos(1),c_pos(2)));
				plot(l_x,l_y,'Color','k');
				xI = xI + 1;
				yI = yI + 1;
			end
		end
	end
	
	c_pos(1) = state.stSubX;
	c_pos(2) = 0;
	for k = 1:length(state.stSubAtoms);
		plot(c_pos(1),c_pos(2),'o');
		[c_x, c_y] = circ(c_pos(1),c_pos(2),radius);
		plot(c_x,c_y,':');
		c_pos(1) = state.stSubX + state.stSubAtoms(k);
		c_pos(2) = c_pos(2);
	end
	
	%refreshdata
	if i == 1
		pause
	else
		if pc ~= -1
			pause(pc);
		else
			pause;
		end
	end

	fprintf('Percent: %.2f Time: %.2f\n',((i/length(y(:,1))) * 100),t(i));

	
	if i ~= length(y(:,1))
		clf(fign);
	end
end


end
