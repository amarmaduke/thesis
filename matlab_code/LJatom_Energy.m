function y = LJatom_Energy( x, epsi, sigma )

	p1 = sigma./x;
	p2 = p1.*p1;
	p4 = p2.*p2;
	p6 = p4.*p2;
	p12 = p6.*p6;
	
	y = epsi*(p12-2*p6);
	%y = epsi*((sigma/x)^12 - 2*(sigma/x)^6);
	%y = epsi*(x-sigma)^2;

end
