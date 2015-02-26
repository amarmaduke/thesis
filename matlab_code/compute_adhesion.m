function [mvForce, mvCount, stForce, stCount] = compute_adhesion( y, state )
%ADHESION Summary of this function goes here
%   Detailed explanation goes here

	mvCount = 0;
	stCount = 0;
	mvForce = 0;
	stForce = 0;
	
	for j = 2:state.N
		adhered = false;
		for i = 1:length(state.mvSubAtoms)
			mvSub_xI = 1;
			mvSub_x = y(mvSub_xI);
			mvSub_y = y(mvSub_xI + state.N);
			
			xps = y(j) - (mvSub_x + state.mvSubAtoms(i));
			yps = y(j+state.N) - mvSub_y;
			d = xps.^2 + yps.^2;
			dist = sqrt(d);
			tolerance = nthroot(2,6)*state.sigma + 1e-6;
			if dist <= tolerance && ~adhered
				mvCount = mvCount + 1;
				adhered = true;
			end
			
			temp_y = yps/dist;
			temp_x = xps/dist;
			LJval = LJatom(dist,state.epsi(3),state.sigma);
			mag = sqrt((LJval*temp_y).^2 + (LJval*temp_x).^2);
			mvForce = mvForce + mag;
		end
	end
	
	for j = 2:state.N
		adhered = false;
		for i = 1:length(state.stSubAtoms)
			xps = y(j) - (state.stSubX + state.stSubAtoms(i));
			yps = y(j+state.N) - 0;
			d = xps.^2 + yps.^2;
			dist = sqrt(d);
			tolerance = nthroot(2,6)*state.sigma + 1e-6;
			if dist <= tolerance && ~adhered
				stCount = stCount + 1;
				adhered = true;
			end
			
			temp_y = yps/dist;
			temp_x = xps/dist;
			LJval = LJatom(dist,state.epsi(2),state.sigma);
			mag = sqrt((LJval*temp_y).^2 + (LJval*temp_x).^2);
			stForce = stForce + mag;
		end
	end
end

