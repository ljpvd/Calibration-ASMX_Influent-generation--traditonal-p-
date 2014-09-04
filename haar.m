% haar wavelet
% Author: Lluis Ma Bosch (6/3/2014) (holalluis@gmail.com)

function [y]=haar(x)

for i=1:length(x)

	if(x(i)>=0 & x(i)<0.5)

		y(i)=1;

	elseif(x(i)>=0.5 & x(i)<1)

		y(i)=-1;

	else

		y(i)=0;

	end
end

end
