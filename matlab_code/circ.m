function [x,y] = circ( x0, y0, r )

	t = -pi:.01:pi;
	x = r*cos(t) + x0;
	y = r*sin(t) + y0;

end

