function[absorb, absorb_lp, t_SW, t_incident, t_EW, spike_ind, comments]=MlappLasTrace(LaserName,SigmaRelations,react,Pretrigger,time, pressureVolts, Pitch, Catch, PitchDark,CatchDark,PxdPsiToVolts,timetoendwall, timefromendwall, LowPassFrequency, OscopeSampleFrequency, comments)
%Calculate absorbance time history and find reflected schlieren spike for LAS laser. Can change which one is used by setting t_SW to the appropriate variable
Ru=8.314;
%Shock Data File
%Pretrigger=str2double(Pretrigger); %[ms]
PretriggerIndex=(Pretrigger/1000)*(OscopeSampleFrequency);
%time = shock_columns(:,1);
Pressure_Pa=pressureVolts .*PxdPsiToVolts*6894.76;% [pa]... %KistlerSpecs(3) = psi/volt
%Pitch = shock_columns(:,3);%Iref-signal (dark offset applied)
PitchAvg=mean(Pitch);
PitchStdev=std(Pitch);
%Catch = shock_columns(:,4);%It-signal (dark offset applied)
%PitchDark = dark_columns(:,1);
PitchDarkAvg = mean(PitchDark);%Iref-vaccuum (dark offset applied)
PitchDarkStdev=std(PitchDark);
%CatchDark = dark_columns(:,2);
CatchDarkAvg = mean(CatchDark);%I0-vacuum (dark offset applied)
CatchDarkStdev=std(CatchDark);

IoBias=CatchDarkStdev; %Claims bias is just std dev of signal
It=Catch.*(PitchDarkAvg/CatchDarkAvg); %It CMR-corrected
absorb = -log((Catch./Pitch).*PitchDarkAvg/CatchDarkAvg);% With Common mode rejection (Dark signal offset applied in LabView)

[~,spike_ind] = max(absorb(1:ceil(length(absorb)/2)));
t_A = time(spike_ind);
absorb_lp = lowpass(absorb,LowPassFrequency,OscopeSampleFrequency);
dAbsdt=gradient(absorb_lp,time);% point-by-point derivative of absorb_lp
absorbBaseline=mean(absorb_lp(1:PretriggerIndex/2));
absorbstd=std(absorb_lp(1:PretriggerIndex/2));
%dAbsdtcontainsImaginary=~isreal(dAbsdt);

if abs(-2*absorbstd)<abs(absorbBaseline)||absorbBaseline>(2*absorbstd) %If baseline absorbance is far from zero, correct it...
    absorb_lp=absorb_lp-absorbBaseline; %Correct to zero so you can still use this to find State 2 & 5 peaks, but:
    disp(strcat('Baseline',LaserName,' absorbance was < 0. Signal might be unreliable'))
    comments=strcat(comments,'Baseline',LaserName,' absorbance was < 0. ');
end %If baseline absorbance is far from zero, don't use absorbance data, but use the schlieren spikes to ID states 2 and 5.

%{
dAbsdtcontainsImaginary=~isreal(dAbsdt);
if dAbsdtcontainsImaginary
    disp('dAbs/dt contains imaginary numbers. The signal was probably bad, and all your LAS-derived data will be written as NA');
    State5Sigma='NA'; State5SigmaBias='NA'; t_SW = 'NA';t_incident='NA';t_EW ='NA';
    State2Abs='NA';State5Abs='NA';State2AbsStdev='NA';PitchDarkAvg='NA';PitchDarkStdev='NA';CatchDarkAvg='NA';CatchDarkStdev='NA';State1Abs='NA';State5AbsStdev='NA';PitchVariation='NA';State2Sigma='NA';State2SigmaBias='NA';State5XBarBias='NA';State5XBarBias='NA';State2XBar='NA';State2XBarBias='NA';CatchVariation='NA';
else
    [pks,locs,w,p] = findpeaks(absorb_lp(1:StopReflectedShockSearch),'MinPeakHeight',2*max(absorb_lp(100:PretriggerIndex/1.5)),'MinPeakDistance',((timetoendwall+timefromendwall)*OscopeSampleFrequency)/2); %pks=peaks, locs=locations, w=width, p=prominence, 'MinPeakHeight' means ignore all peaks below 1500.
    State1End=locs(1)-int64(((timetoendwall+timefromendwall)*OscopeSampleFrequency)/1);%index of State 1 start
    State2Start=locs(1)+int64(((timetoendwall+timefromendwall)*OscopeSampleFrequency)/10); %index of State 2 start
    State2End=locs(find(pks==max(pks)))-int64(((timetoendwall+timefromendwall)*OscopeSampleFrequency)/3); %index of State 2 end
    Catch1= mean(shock_columns(4:State1End,3));

    %Catch signal variation suggests that either the laser is unstable, or there is significant background IR emitters.
    CatchVariation=CatchDarkAvg/Catch1;%Catch Detector Variation
    PitchVariation=PitchDarkAvg/PitchAvg;%Pitch Detector Variation
    if CatchVariation<0.95||CatchVariation>1.05 %If there's significant variation in the Catch detector....
        if PitchVariation<0.95||PitchVariation>1.05%...And significant variation in the Pitch detector...
            if PitchVariation<(0.99*CatchVariation)||PitchVariation>(1.01*CatchVariation) %...and the variation in the two is similar...
                disp(strcat(LaserName,' PITCH Vac = ', num2str(PitchVariation),'x PITCH Sample, and ',LaserName,' CATCH Vac = ', num2str(CatchVariation),'x ',LaserName,' State 1 Sample. Its likely the laser was unstable.'));
                comments=strcat(comments,strcat(LaserName,' PITCH Vac = ', num2str(PitchVariation),'x PITCH Sample, and ',LaserName,' CATCH Vac = ', num2str(CatchVariation),'x',LaserName,' State 1 Sample. Its likely the laser was unstable.'));
            else %But if the variation in the two detectors is dissimilar...
                disp(strcat(LaserName,' PITCH Vac = ', num2str(PitchVariation),'x PITCH Sample, and ',LaserName,' CATCH Vac = ',num2str(CatchVariation),'x ',LaserName,' State 1 Sample. Its likely the laser was unstable AND there were background IR emitters.'));
                comments=strcat(comments,strcat(LaserName,' PITCH Vac = ', num2str(PitchVariation),'x PITCH Sample, and ',LaserName,' CATCH Vac = ', num2str(CatchVariation),'x ',LaserName,' State 1 Sample. Its likely the laser was unstable AND there were background IR emitters.'));
            end
        else %If there's variation in the Catch detector, but not the Pitch detector...
            strcat(LaserName,' PITCH Vac = ', num2str(PitchVariation),'x PITCH Sample, But ',LaserName,' CATCH Vac = ',num2str(CatchVariation),'x ',LaserName,' State 1 Sample. Its likely the laser was stable but there were background IR emitters.')
            comments=strcat(comments,strcat(LaserName,' PITCH Vac = ', num2str(PitchVariation),'x PITCH Sample, and ',LaserName,' CATCH Vac = ', num2str(CatchVariation),'x ',LaserName,' State 1 Sample. Its likely the laser was stable but there were background IR emitters.'));
        end
    end

    State5start=locs(find(pks==max(pks)))+int64(((timetoendwall+timefromendwall)*OscopeSampleFrequency)/5); %index of State 5 start

    if abs(-2*absorbstd)<abs(absorbBaseline)||absorbBaseline>(2*absorbstd) %If baseline absorbance is far from zero, don't use the absorbance data.
        State1Abs='NA'; State2Abs='NA'; State2AbsStdev='NA'; State2XBar='NA'; State5Abs='NA'; State5AbsStdev='NA'; State2Sigma='NA'; State2SigmaBias='NA';  State5Sigma='NA'; State5SigmaBias='NA'; State5XBarBias ='NA'; State5XBarBias='NA';
        disp('Baseline LAS absorbance 2*std away from 0. LAS Abs NA.')
        comments=strcat(comments,'BL Abs_LAS > 2*std away from 0. ');State5end=1;
    else
        State1Abs=mean(absorb(1:State1End));
        State2Abs=mean(absorb(State2Start:State2End));
        State2AbsStdev=std(absorb(State2Start:State2End));
        State2It=mean(It(State2Start:State2End));
        State2ItBias=std(It(State2Start:State2End));%Claims bias is just std dev of signal
        if strcmpi(react,'y') %Reacting Tests:  state 5 ends when ignition begins.
            %    if strcmpi(UserInput{6},'y')
            %         State5end=(index_intercept_OH_EW)-500; % If OH endwall was used, use that, cuz it should happen earliest
            %         disp('OH Endwall tign used for CO2 State 5 end')
            %         comments=strcat(comments,'OH EW tign used for CO2 State 5 end');
            %    else
            %        if strcmpi(UserInput{7},'y') %if OH endwall wasn't used, but OH sidewall was,
            %           State5end=(index_intercept_OH_SW)-500; % Use OH sidewall, cuz it's more reliable than pressure
            %           disp('OH Sidewall tign used for CO2 State 5 end')
            %           comments=strcat(comments,'OH SW tign used for CO2 State 5 end');
            %        else
            State5end=(t_ignPIndex)-500; % If you have not other choice, use P
            disp('PXD tign used for',LaserName,' State 5 end')
            comments=strcat(comments,'PXD tign used for ',LaserName,' State 5 end');
            %        end
            %    end
            State5Abs=mean(absorb(State5start:State5end));
            State5AbsStdev=std(absorb(State5start:State5end));
            State5It=mean(It(State5start:State5end));
            State5ItBias=std(It(State5start:State5end));%Claims bias is just std dev of signal
            State2Sigma='NA';State2SigmaBias='NA';State5Sigma='NA'; State5SigmaBias='NA';
        else %Non-reacting (Absorbing) tests
            %Rarefaction_index = find(dPdt < 0);
            %State5end=Rarefaction_index(1);  %Rarefaction identification still not working right
            State5end=State5start+int64(((timetoendwall+timefromendwall)*OscopeSampleFrequency)*4);
            disp('State 5 start State2 Lengthx4  used to determine CO2 Laser State 5 end.  Still need to switch to using Rarefaction wave arrival')
            comments=strcat(comments,'State 5 start State2 Lengthx4 used to determine CO2 State 5 end');
            State5Abs=mean(absorb(State5start:State5end));
            State5AbsStdev=std(absorb(State5start:State5end));
            State5It=mean(It(State5start:State5end));
            State5ItBias=std(It(State5start:State5end));%Claims bias is just std dev of signal
            State2Sigma=(State2Abs*Ru*Tout(1))/(XtestGas(1)*Pout(1)*TubeDiameter);
            State2DelSigmaDelX=(-State2Abs*Ru*Tout(1))/(XtestGas(1)^2*TubeDiameter*Pout(1)); %influence coefficient - X_C2H4
            State2DelSigmaDelP=(-State2Abs*Ru*Tout(1))/(XtestGas(1)*TubeDiameter*Pout(1)^2); %influence coefficient - P2
            State2DelSigmaDelL=(-State2Abs*Ru*Tout(1))/(XtestGas(1)*TubeDiameter^2*Pout(1)); %influence coefficient - L
            State2DelSigmaDelT=(-State2Abs*Ru)/(XtestGas(1)*TubeDiameter*Pout(1)); %influence coefficient - T2
            State2DelSigmaDelIt=(-Ru*Tout(1))/(XtestGas(1)*TubeDiameter*Pout(1)*State2It); %influence coefficient - It
            State2DelSigmaDelI0=(-Ru*Tout(1))/(XtestGas(1)*TubeDiameter*Pout(1)*CatchDarkAvg); %influence coefficient - I0
            State2SigmaBias=sqrt((State2DelSigmaDelX*XtestGas_bias(1))^2+(State2DelSigmaDelP*P2_bias)^2+(State2DelSigmaDelL*Bias_TubeDiameter)^2+(State2DelSigmaDelT*T2_bias)^2+(State2DelSigmaDelIt*State2ItBias)^2+(State2DelSigmaDelI0*IoBias)^2);  %[m^2/mol]
            State5Sigma=(State5Abs*Ru*Tout(2))/(XtestGas(1)*Pout(2)*TubeDiameter);
            State5DelSigmaDelX=(-State5Abs*Ru*Tout(2))/(XtestGas(1)^2*TubeDiameter*Pout(2)); %influence coefficient - X_C2H4
            State5DelSigmaDelP=(-State5Abs*Ru*Tout(2))/(XtestGas(1)*TubeDiameter*Pout(2)^2); %influence coefficient - P2
            State5DelSigmaDelL=(-State5Abs*Ru*Tout(2))/(XtestGas(1)*TubeDiameter^2*Pout(2)); %influence coefficient - L
            State5DelSigmaDelT=(-State5Abs*Ru)/(XtestGas(1)*TubeDiameter*Pout(2)); %influence coefficient - T2
            State5DelSigmaDelIt=(-Ru*Tout(2))/(XtestGas(1)*TubeDiameter*Pout(2)*State5It); %influence coefficient - It
            State5DelSigmaDelI0=(-Ru*Tout(2))/(XtestGas(1)*TubeDiameter*Pout(2)*CatchDarkAvg); %influence coefficient - I0
            State5SigmaBias=sqrt((State5DelSigmaDelX*XtestGas_bias(1))^2+(State5DelSigmaDelP*P5_bias)^2+(State5DelSigmaDelL*Bias_TubeDiameter)^2+(State5DelSigmaDelT*T5_bias)^2+(State2DelSigmaDelIt*State5ItBias)^2+(State5DelSigmaDelI0*IoBias)^2);  %[m^2/mol]
        end
    end
    [State2Sigma, State2SigmaBias]=SigmaCorrelations(Tout(1),SigmaRelations);
    [State5Sigma, State5SigmaBias]=SigmaCorrelations(Tout(2),SigmaRelations);
    [StateiSigma, StateiSigmaBias]=SigmaCorrelations(Ti,SigmaRelations);

    if strcmpi(SigmaRelations,'None') %If there's no relation to allow you to determine sigma
        State2XBar='NA';State2XBarError='NA'; State2XBarBias ='NA'; State5XBar='NA';State5XBarError='NA'; State5XBarBias ='NA';StateiXBar='NA';
    else %Otherwise, if there is a relation available, use it, and calculate/plot species mole fraction.
        State2XBar=(State2Abs*Ru*Tout(1))/(State2Sigma*Pout(1)*TubeDiameter);  %Single point, average value
        State2XBarError=1-State2XBar/XtestGas(1); %percent error compared to manometric measurement
        State2XBarBias ='NA';  %not calculated YET.  working on it.
        State5XBar=(State5Abs*Ru*Tout(2))/(State5Sigma*Pout(2)*TubeDiameter);
        State5XBarError=1-State5XBar/XtestGas(1); %percent error compared to manometric measurement
        State5XBarBias ='NA';  %not calculated YET.  working on it.
        StateiXBar=(absorb_lp.*Ru.*Ti)./(StateiSigma.*Pressure_Pa.*TubeDiameter);  %calculatated X at every point

        figureIndex=figureIndex+1;
        figure(figureIndex)
        plot(StateiXBar,'DisplayName',strcat('X ',LaserName,' raw'))
        hold on
        yline(State2XBar,'green--','DisplayName',strcat('State2 X ',LaserName))
        yline(State5XBar,'r-','DisplayName',strcat('State5 X ',LaserName))
        yline(XtestGas(1),'DisplayName',strcat('Manometric X ',LaserName))
        legend
        hold off
    end

 
    figureIndex=figureIndex+1;
    figure(figureIndex)
    hold on
    plot(time, Pitch./mean(Pitch(1:PretriggerIndex/2)),'r-','DisplayName',strcat(LaserName,' Ref Norm'))
    plot(time, Catch./Catch1,'b-','DisplayName',strcat(LaserName,' It Norm'))
    %plot(time,CatchDark./Catch1,'b-.','DisplayName',strcat(LaserName,'I0 Norm')) %This one in particular makes Matlab choke.
    plot(time, absorb,'k-','DisplayName',strcat(LaserName,' absorbance'))
    plot([time(State1End) time(State1End)], [0 1],'green--','DisplayName','State1End')
    plot([time(State2Start) time(State2Start)], [0 1],'k--','DisplayName','State2Start')
    plot([time(State2End) time(State2End)], [0 1],'k-','DisplayName','State2End')
    plot([time(State5start) time(State5start)], [0 1],'k:','DisplayName','State5start')
    plot([time(State5end) time(State5end)], [0 1],'k-.','DisplayName','State5End')
    xlabel('time [ms]')
    ylabel('Signal [-]')
    legend
    xlim([1.5 5])
    hold off

    figureIndex=figureIndex+1;
    figure(figureIndex)
    plot(absorb,'DisplayName',strcat('1)absorb ',LaserName,' raw'))
    hold on
    plot(absorb_lp,'DisplayName',strcat('2)absorb ',LaserName,'  lp'))
    yline(absorbBaseline,'DisplayName','Real Baseline Abs')
    plot(dAbsdt./100,'DisplayName',strcat('4)d',LaserName,'dt./1000'))
    plot(locs,pks,'o','MarkerSize',12,'DisplayName','2)Inc & Ref Shock pks')
    legend
    hold off
%}

t_SW = t_A; %[ms] time of REFLECTED? shock arrival at sidewall.
t_incident ='NA';% time(locs(1));%[ms] time of INCIDENT shock arrival at ENDWALL?
t_EW ='NA';% t_incident + (timetoendwall/1000);%[ms] time of incident shock arrival at end wall

%AbsParameters=[State2Abs,State5Abs,State2AbsStdev,PitchDarkAvg,PitchDarkStdev,CatchDarkAvg,CatchDarkStdev,PitchAvg, PitchStdev,State1Abs,State5AbsStdev,PitchVariation, State2Sigma, State2SigmaBias, State5XBarBias, State5XBarBias,State2XBar, State2XBarBias, CatchVariation];
end