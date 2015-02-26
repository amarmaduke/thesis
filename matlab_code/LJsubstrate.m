function y = LJsubstrate( x, epsi, sigma, density )

	p1 = sigma/x;
	p2 = p1*p1;
	p4 = p2*p2;
	p5 = p4*p1;
	p11 = p5*p5*p1;
	
	y = -(pi*epsi*density)*(2*p11-4*p5);
%    y=100*(x-1);
end

