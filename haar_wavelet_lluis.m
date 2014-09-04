function yapp=haar_wavelet_lluis(y,T)
	
	% y: noisy vector
	% T: distancia de salt

	x=1:length(y);

	%nova figura
	figure(),grid,hold on

		%original data
		plot(x,y,'LineWidth',2)

		%filtered data
		yapp=[];

		for i=x(1):T:x(end)
			inici=i;
			final=inici+T;
			p_mig=(inici+final)/2;
			% arrays logics
			interval___ = (x>=inici&x<=final);
			interval_Lo = (x>=inici&x<=p_mig);
			interval_Up = (x>=p_mig&x<=final);
			C___ = mean(y(interval___));
			C_Lo = mean(y(interval_Lo));
			C_Up = mean(y(interval_Up));
			scaling=(C___-C_Lo);
			yapp(interval___)=scaling*-haar((x(interval___)-inici)/T)+C___;
		end

		plot(x,yapp,'r','LineWidth',2);
		legend('dades originals',['haar wavelet noise reduction, T=',num2str(T)])

end
