function y = LJsubstrate_Energy( x, epsi, sigma, density )

	p1 = sigma/x;
	p2 = p1.*p1;
	p4 = p2.*p2;
	p8 = p4.*p4;
	p10 = p8.*p2;
	
	y = (pi*epsi*sigma*density^2)*(.2*p10-p4);
%    y=50*(x-1)^2;
end
