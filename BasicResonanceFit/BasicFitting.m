function  [v_mag,oml,om0,mag,pha,Gain] = BasicFitting(pathToFile,offs,ang,ph)

    load(pathToFile, 'f','mag','pha','Gain','SignalGeneratorPumpP');

    freqs_exper=f-f(length(f)/2+0.5);
    om0=f(length(f)/2+0.5);
    pows_exper=SignalGeneratorPumpP;
    pows_exper(1)=[];
    Gain(:,1)=[];
    magn=db2mag((mag(:,:)-mag(:,1)));
    magn(:,1)=[];
    mag_init=db2mag(mag(:,2).'-mag(:,1).')-0.1;
    pha_init=pha(:,2).'-pha(:,1).'+ph;
    pha(:,:)=pha(:,:)-pha(:,1)+ph;
    pha(:,1)=[];
    om0=2*pi*om0;
    oml=linspace(freqs_exper(1),freqs_exper(end),length(freqs_exper))+offs;%+0.0683e6;
    om=oml*2*pi;
   
    fit_func0 = @(v,oml)(1-sqrt(1)*(exp(1i*ang*pi)*(v(1))./(-1i*2*pi*(oml)+(v(1)+v(2))/2)));

    opts = optimoptions('lsqcurvefit','Algorithm','trust-region-reflective');
    opts.OptimalityTolerance=1e-12;
    opts.MaxIterations=1600;
    opts.FunctionTolerance=1e-12;
    opts.StepTolerance=1e-12;
%     ph=0;
    x0 = [1e7;1e7]; 
    [v_mag,resnorm] = lsqcurvefit(fit_func0,x0,oml,(((mag_init).*exp(1i.*(pha_init)))),[],[],opts)

    figure (52)
    plot(freqs_exper+offs,[abs(1+(exp(1i*-ang*pi).*(-1+(mag_init.'.*exp(1i.*(pha_init.')))))),angle(1+(exp(1i*-ang*pi).*(-1+(mag_init.'.*exp(1i.*(pha_init.'))))))]);
    hold on
    plot(oml,[(abs((1+exp(1i*0*pi)*v_mag(1)./(1i*oml*2*pi-(v_mag(1)+(v_mag(2)))/2))));angle((1+exp(1i*0*pi)*v_mag(1)./(1i*oml*2*pi-(v_mag(1)+(v_mag(2)))/2)))],'.')
    hold off
    xlabel('Frequency (MHz)');
    ylabel('Normalized magnitude (a.u.) and phase (rad)');
    legend('Mag Measurement','Ph Measurement','Mag Fitting','Ph Fitting')
end

