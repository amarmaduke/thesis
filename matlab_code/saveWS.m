path = './';
clck = clock;
name = strcat(num2str(clck(1)),'-',num2str(clck(2)),'-', ...
			num2str(clck(3)),',',num2str(clck(4)),':',num2str(clck(5)));
str = strcat(path,name,'.mat');
save(str,'t','y','state');