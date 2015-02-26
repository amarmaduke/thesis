function [ j, i ] = getAtom( x, n )
		temp_s = n(1);
		j = 1;
		for k = 2:length(n)
			if x <= temp_s
				temp_s = temp_s - n(k-1);
				break;
			end
			j = k;
			temp_s = temp_s + n(k);
		end
		if j == length(n)
			temp_s = temp_s - n(length(n));
		end
		i = x - temp_s;
end

