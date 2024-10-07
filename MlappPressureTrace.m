function[PlotData, t_ignPIndex, IDT_P, IDT_P_bias]=MlappPressureTrace(Reacting, Pretrigger, timetoendwall, timefromendwall, dark_columns,time, Kistler_Voltage,KistlerSpecs,Uis,OscopeSampleFrequency)

%Pretrigger=str2double(UserInput{11}); %[ms]
PretriggerIndex=(Pretrigger/1000)*(OscopeSampleFrequency);
pressurenoise=dark_columns(:,1);
Pressure_atm_raw=Kistler_Voltage.*(KistlerSpecs(3)/14.7);% [atm]... %KistlerSpecs(3) = psi/volt                       
disp('Temperature for CO2 absorbance calculations based on constant volume assumption. I **Think** it should be adiabatic, where PV^n (n=gamma) is constant, not volume..')
PXD_ShockPass_Frequency=Uis/KistlerSpecs(2); %KistlerSpecs(2) = PXD diameter[m]
PXD_fft_Frequency=55000; %Prove to yourself this is a common noise frequency by uncommenting the next two lines of code.
%fftP=fft(pressurenoise);
%plot(OscopeSampleFrequency/length(pressurenoise)*(-length(pressurenoise)/2:length(pressurenoise)/2-1),abs(fftshift(fftP)),"LineWidth",3)
PXD_LP_Frequency=min([PXD_ShockPass_Frequency KistlerSpecs(1) PXD_fft_Frequency]);
pressureSG= sgolayfilt(Kistler_Voltage,2,51);
BWorder = 2; % Order of the Butterworth filter
[b, a] = butter(BWorder, PXD_LP_Frequency/ (OscopeSampleFrequency/ 2), 'low');% Design Butterworth filter
pressureBW = filtfilt(b, a, pressureSG);% Apply filter using filtfilt to achieve zero-phase filtering
PressureMerge = sqrt(abs(pressureBW .* pressureSG));
pressureLP = PressureMerge;% lowpass(Kistler_Voltage,PXD_LP_Frequency,OscopeSampleFrequency); %y = lowpass(x,fpass,fs)
Pressure_atm=pressureLP*KistlerSpecs(3)*0.068046;
dPdt=gradient(pressureLP,time);% point-by-point derivative of absorbCO2
%dPdtSG= sgolayfilt(dPdt,2,51);
%dPdtBW = filtfilt(b, a, dPdtSG);
%dpdtMerge =sqrt(abs(dPdtBW .* dPdtSG));
%dpdt=dPdtSG;
PressureWindowCenter=int64(PretriggerIndex); %[index]center of "findpeaks" window.
ReflectedToIncidentBackwardsIndex=int64((timetoendwall+timefromendwall)*(OscopeSampleFrequency));
%StartIncidentShockSearch=int64(PressureWindowCenter-(ReflectedToIncidentBackwardsIndex*2));
%StopReflectedShockSearch=int64(PressureWindowCenter+3000);%Change 3000 to some a shock-specific number.
%{
%Find incident and reflected shock peaks
%[Ppks,Plocs,Pw,Pp] = findpeaks(dPdt(StartIncidentShockSearch:StopReflectedShockSearch),'MinPeakHeight',0.3*max(dPdt));%,'Npeaks',2,'MinPeakDistance',length(StartIncidentShockSearch:StopReflectedShockSearch)/2); %pks=peaks, locs=locations, w=width, p=prominence, 'MinPeakHeight' means ignore all peaks below 1500.
[Ppks,Plocs,Pw,Pp] = findpeaks(dPdt(StartIncidentShockSearch:StopReflectedShockSearch),'MinPeakProminence',100,'Npeaks',2,'MinPeakDistance',((timetoendwall+timefromendwall)*OscopeSampleFrequency)/2); %pks=peaks, locs=locations, w=width, p=prominence, 'MinPeakHeight' means ignore all peaks below 1500.
Plocs=[Plocs(1)+StartIncidentShockSearch Plocs(2)+StartIncidentShockSearch];%redefine Plocs in terms of master index.
Ppks=[Ppks(1) Ppks(2)]; %make it same length as Plocs, in case there were more than 2 peaks.
PState1End=Plocs(1)-int64(((timetoendwall+timefromendwall)*OscopeSampleFrequency)/1);%index of State 1 start 
PState2Start=Plocs(1)+int64(((timetoendwall+timefromendwall)*OscopeSampleFrequency)/10); %index of State 2 start 
PState2End=Plocs(2)-int64(((timetoendwall+timefromendwall)*OscopeSampleFrequency)/3); %index of State 2 end
PState5start=Plocs(2)+int64(((timetoendwall+timefromendwall)*OscopeSampleFrequency)/5); %index of State 5 start
t_A_PXD = time(Plocs(2)); %but this won't work yet cuz your PXD filtering sucks.
%}

%Find ignition from P
if strcmpi(Reacting,'y') %Reacting Tests
    [P_max,P_max_ind] = max(pressureLP);
    [dPdtMax,dPdtMaxIndex] = max(dPdt(P_max_ind-3000:P_max_ind)); %Need to improve this too....point of maximum slope in the 3000 data points leading up to max P.
    t_ignPIndex=P_max_ind-3000+dPdtMaxIndex;
    t_ignP=time(t_ignPIndex);
    t_SW_P = 'Na';% time(Plocs(2));%[ms] %use max
    IDT_P= 'Na';% (t_ignP-t_SW_P)*1000; %[us]  
    IDT_P_bias='NA';  %not calculated YET.  working on it.
    PState5end=(t_ignPIndex); % state 5 ends when ignition begins. 
else %Non-reacting (Absorbing) tests
    Rarefaction_index = find(dPdt < 0);
    %State5end=Rarefaction_index(1);  %Rarefaction identification still not working right
    PState5end=PState5start+int64(((timetoendwall+timefromendwall)*OscopeSampleFrequency)*4);
    disp('State 5 start + State2 Lengthx4 used to determine Kistler State 5 end.  Still need to switch to using Rarefaction wave arrival')
   t_ignPIndex='NA';t_SW_P='NA'; IDT_P='NA';IDT_P_bias='NA';
end

PlotData=[time, Pressure_atm, Pressure_atm_raw, pressureLP];
figure
plot(time,Pressure_atm,'DisplayName',"Pressure [atm]")
hold on
plot(time,Pressure_atm_raw,'DisplayName',"Pressure raw [atm]")
plot(time,pressureLP,'DisplayName',"LP Pressure [atm]")
legend
ylabel="Pressure [atm]";
xlabel('time [ms]')

%{
figureIndex=figureIndex+1;
figure(figureIndex)
hold on
yline(Pout(2)/101325,'DisplayName','FROSH P5')
plot(time,Pressure_atm,'DisplayName','PXD')
if strcmpi(Reacting,'y')
plot([t_ignP t_ignP], [0 max(Kistler_Voltage)*20/14.7],'b--','DisplayName','t Ign P')
else
end
%{
plot([time(PState1End) time(PState1End)], [0 max(Kistler_Voltage)*20/14.7],'green--','DisplayName','State1End')
plot([time(PState2Start) time(PState2Start)], [0 max(Kistler_Voltage)*20/14.7],'k--','DisplayName','State2Start')
plot([time(PState2End) time(PState2End)], [0 max(Kistler_Voltage)*20/14.7],'k-','DisplayName','State2End')
plot([time(PState5start) time(PState5start)], [0 max(Kistler_Voltage)*20/14.7],'k:','DisplayName','State5start')
plot([time(PState5end) time(PState5end)], [0 max(Kistler_Voltage)*20/14.7],'k-.','DisplayName','State5End')
%}
xlabel('time [ms]')
ylabel('Pressure [atm]')
legend
xlim([1.5 5])
ylim([0 max(Kistler_Voltage)*20/14.7])
hold off

figureIndex=figureIndex+1;
figure(figureIndex)
hold on
plot(Kistler_Voltage,'DisplayName','1)pressure raw')
plot(pressureLP,'DisplayName','2)pressure lp')
plot(dPdt./100,'DisplayName','4)dPdt./1000')
%plot(Plocs,Ppks./100,'o','MarkerSize',12,'DisplayName','2)Inc & Ref Shock pks')
xlim([StartIncidentShockSearch StopReflectedShockSearch])
legend
hold off
%}