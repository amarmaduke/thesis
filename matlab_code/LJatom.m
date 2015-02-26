function y = LJatom( x, epsi, sigma )

	p1 = sigma/x;
	p2 = p1*p1;
	p4 = p2*p2;
	p7 = p4*p2*p1;
	p8 = p7*p1;
	p13 = p8*p4*p1;
	
	y = -(12*epsi/sigma)*(p13-p7);
	%y = -(12*epsi/sigma)*((sigma/x)^13 - (sigma/x)^7);
	%y = 2*epsi*(x-sigma);
end
