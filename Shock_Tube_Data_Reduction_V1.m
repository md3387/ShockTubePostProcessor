%Shock Tube Data Reduction - Immediate DAQ Computer Postprocessing version
%Original version:  MSK Fall 2023

%Function Dependencies:
%frosh_bias_V2.m
%CHON_MW.m
%ParseElementString.m

%Revision History: 
%MDH 11/20/2023 - Added Uis Plot
%MDH 11/21/2023 - Added user dialog & handling of 1-4 mixture components
%MDH 11/22/2023 - Added lots of plotting capabilities.
%MDH 11/28/2023 - Added dialog to use/ignore test section diagnostic
%MDH 11/28/2023 - Added IDT from Pressure
%MDH 12/1/2023 - Boatload of Changes. Renamed to _V1

%Future upgrades
%1.a Change user dialog to look more like a table of values so users don't have to format their inputs
%1.b Can you make dialog defaults just equal the last entry?
%1.d update frosh_bias.m to handle correllated, n-component mixtures - might require 1.e
%1.e Use Monte Carlo analysis to determine contribution of correlated mixture uncertainties
%3.a output State 2 variables (aka absorbance)
%3.c Xinitial is always known. Use this to plot X time-history?
%3.d use all previous D/F data from the current day to plot I0 over time.
%4.b allow user to remove "bad" counter-timer outputs from Uis determination
%5.a figure out how to adjust IDT calculation parameters to eliminate matlab errors
%6.a Make each diagnostic a function, and let user select the functions they want to call, based on the diagnostics

close all
clearvars
%clc
%Define sampling and filtering parameters
OscopeSampleFrequency=10000000; %[Hz]
PhotodetectorFrequency=2600; %[Hz] sidewall
PMTFrequency=80000000; %[Hz] end wall
PXD_Amp_Frequency=100000;%[Hz] Kistler amplifier sample frequency
Kistler_psi_per_volt=20;%[psi/volt]

Ru=8.314; %J/mol*K

%Define file naming parameters
file_base = 'C:\Users\Shock Tube DAQ\Documents\Shock Tube Experiment Data\'; %Location of folder containing all data files
UserPrompt={'Date','Run','Driven Gas Components', 'Driven Gas Raw Baratron Rdgs[Torr]', 'Driver','OH Endwall? y/n','OH Sidewall? y/n','CO2 Laser? y/n','HeNe Laser? y/n','Reacting? y/n','Pretrigger [ms]'};
%Promptdefaults = {'202311dd','001','C2H4 O2 N2','32.72 131 500', 'He','y','y','y','y','y'};
Promptdefaults = {'202310dd','001','C2H4 O2 N2','32.72 131 500', 'He','y','y','y','y','y','2'};
UserInput = inputdlg(UserPrompt,'Input',1,Promptdefaults);
date = UserInput{1}; %Experiment date
run_no = UserInput{2}; %Experiment number, ###
testGasSpec = textscan(UserInput{3}, '%s'); %parse the string to produce a cell array with three values
testGasSpec =testGasSpec{1}'; %pull values out of cell to produce a linear array that frosh can use
Pmix = textscan(UserInput{4}, '%f');%parse the string to produce a cell array with three values
Pmix = [0 Pmix{1}'];%pull values out of cell to produce a linear array, and add the vacuum pressure (0 Torr) to the front of the array
Pretrigger=str2double(UserInput{11}); %[ms]
PretriggerIndex=(Pretrigger/1000)*(OscopeSampleFrequency);

%Partial Pressure Bias Calclation
for i=1:length(Pmix)
    if Pmix(i)<100 %If you used the 100 Torr Baratron...
        PmixBias(i)=0.01;%[Torr]
    else          %If you used the 10000 Torr Baratron...
        PmixBias(i)=1;%[Torr]
    end
end
for i=1:(length(Pmix)-1)%there is one less mole fraction than pressure reading.
        XtestGas(i)=(Pmix(i+1)-Pmix(i))/Pmix(length(Pmix));
        XtestGas_bias(i)=sqrt(((-1/Pmix(length(Pmix))*PmixBias(i))^2)+((1/Pmix(length(Pmix))*PmixBias(i+1))^2)+((((Pmix(i+1)-Pmix(i))/((Pmix(length(Pmix)))^2))*PmixBias(i+1))^2)); %Used to be: XtestGas_bias = XtestGas(1)*0.00707; 
        MW(i)=CHON_MW(char(testGasSpec(i)));
end
MW_mix = sum(XtestGas.*MW);
comments='';
%%


%Tube Physical Dimensions
DM1_DM2=12.00*0.0254;% [m] Dead-mass ports center-to-center spacing (measured. See Table9 (P.34) Ref.1)
DM2_DM3=11.993*0.0254;% [m] Dead-mass ports center-to-center spacing (measured. See Table9 (P.34) Ref.1)
DM3_DM4=11.996*0.0254;% [m] Dead-mass ports center-to-center spacing (measured. See Table9 (P.34) Ref.1)
DM4_DM5=11.994*0.0254;% [m] Counter0 spacing (PCB 1-2) Dead-mass ports center-to-center spacing (measured. See Table9 (P.34) Ref.1)
DM5_DM6=11.991*0.0254;% [m] Counter1 spacing (PCB 2-3) Dead-mass ports center-to-center spacing (measured. See Table9 (P.34) Ref.1)
DM6_DM7=11.992*0.0254;% [m] Counter2 spacing (PCB 3-4) Dead-mass ports center-to-center spacing (measured. See Table9 (P.34) Ref.1)
NumberofCounters=5; % [-] number of separate incident velocity measurements, usually completed by the counter/timers.
Ctr0Dist=DM3_DM4; % [m] %use in solution loop. If ctr locations change, change it here only, and soln loop will automatically reflect the change.
Ctr1Dist=DM4_DM5; % [m] %use in solution loop. If ctr locations change, change it here only, and soln loop will automatically reflect the change.
Ctr2Dist=DM5_DM6; % [m] %use in solution loop. If ctr locations change, change it here only, and soln loop will automatically reflect the change.
Ctr3Dist=DM6_DM7; % [m] %use in solution loop. If ctr locations change, change it here only, and soln loop will automatically reflect the change.
EndwalltoPCB1=61.9*0.0254;% [m] Distance from PCB to endwall (specified. Ref.2 CST-1005, CST-1006, CST-1007)
EndwalltoPCB2=EndwalltoPCB1-DM3_DM4;% [m] Distance from PCB to endwall
EndwalltoPCB3=EndwalltoPCB2-DM4_DM5;% [m] Distance from PCB to endwall
EndwalltoPCB4=EndwalltoPCB3-DM5_DM6;% [m] Distance from PCB to endwall
EndwalltoPCB5=EndwalltoPCB4-DM6_DM7;% [m] Distance from PCB to endwall
EndwalltoTestSection=0.6*0.0254;% Distance from Test Section centerline to endwall (Specified. Ref.2 CST-1006, CST-1007)
EndwalltoMidpoint_DM3_DM4=EndwalltoPCB1-(DM3_DM4/2);% [m] Distance from port spacing midpoint to endwall
EndwalltoMidpoint_DM4_DM5=EndwalltoPCB2-(DM4_DM5/2);% [m] Distance from port spacing midpoint to endwall
EndwalltoMidpoint_DM5_DM6=EndwalltoPCB3-(DM5_DM6/2);% [m] Distance from port spacing midpoint to endwall
EndwalltoMidpoint_DM6_DM7=EndwalltoPCB4-(DM6_DM7/2);% [m] Distance from port spacing midpoint to endwall
EndwalltoMidpoints=[EndwalltoMidpoint_DM3_DM4 EndwalltoMidpoint_DM4_DM5 EndwalltoMidpoint_DM5_DM6 EndwalltoMidpoint_DM6_DM7];
WindowDiam=0.0127;%[meters] = 0.5in -Diameter of test section windows
TubeDiameter=0.1016; %[m]
Bias_TubeDiameter=0.0000127; %[m]
%Load Counter/Timer Data
CT_file = strcat(file_base,num2str(date),'_',num2str(run_no,'%03d'),'_','TimerData.csv');
CT_data_full = readmatrix(CT_file,'NumHeaderLines',1);
CT_data = CT_data_full(:,end);

%Extract Counter/Timer data and calculate incident shock velocity and its uncertainty
U_CT0 = DM3_DM4/CT_data(1);
U_CT1 = DM4_DM5/CT_data(2);
U_CT2 = DM5_DM6/CT_data(3);
U_CT3 = DM6_DM7/CT_data(4);
U=[U_CT0 U_CT1 U_CT2 U_CT3];
weight=[1 1 1 1];

if any(isnan(U))||any(find(U<500))||any(find(U>5000)) %if individual CT measurement returns NAN, U<500m/s, or U>5000m/s, disregard it. 
    if any(isnan(U))
    U_NAN_index=find(isnan(U)); %Index of crappy shock velocity
    elseif any(find(U<500))
    U_NAN_index=find(U<500); %Index of crappy shock velocity
    else
    U_NAN_index=find(U>5000); %Index of crappy shock velocity
    end
    U(U_NAN_index) = [];
    EndwalltoMidpoints(U_NAN_index) = [];
    weight(U_NAN_index) = [];
    U_NAN_LabviewIndex= U_NAN_index-1;% Correct for stupid Labview zero-indexing.
    disp(strcat('C/T',num2str(U_NAN_LabviewIndex), ' was U<500, U>5000, or NAN.  It will be removed from the U_IS calculation, and your uncertainty will increase accordingly.'))
    comments=strcat(comments,'C/T',num2str(U_NAN_LabviewIndex), ' Hi/Lo/NAN. Uis_bias T5_bias P5_bias likely high ');
else
    %do nothing.
end

    %Regular, linear fit (1st degree polynomial)
U_fit = fit(EndwalltoMidpoints',U','poly1');
Uis = U_fit.p2;
timetoendwall=EndwalltoTestSection/Uis; %[s]use along with timefromendwall for estimating expected time between incident and reflected shock phenomena
Uis_CI = confint(U_fit);
Uis_bias = (Uis_CI(2,2)-Uis_CI(1,2))/2;
%Weighted nonlinear regression model, 1st degree polynomial, Weights all =1, so should give same result as above, but easier to plot Confidence intervals
[model, S]=polyfit(EndwalltoMidpoints,U,1);
start=[model(1);model(2)];
modelstring=@(b,EndwalltoMidpoints) b(1).*EndwalltoMidpoints+b(2);
U_wnlm=fitnlm(EndwalltoMidpoints,U,modelstring,start,'Weight',weight);
xx=linspace(0,2)';
[Upred,UpredCI]=predict(U_wnlm,xx,'Simultaneous',true);

figure(2)
hold on
plot(EndwalltoMidpoint_DM3_DM4,U_CT0,'r O')
plot(EndwalltoMidpoint_DM4_DM5,U_CT1,'r O')
plot(EndwalltoMidpoint_DM5_DM6,U_CT2,'r O')
plot(EndwalltoMidpoint_DM6_DM7,U_CT3,'r O')
plot(U_fit)
%%line(xx,predict(U_wnlm,xx),'color','b','linestyle', '--')
plot(xx,Upred,'b--',xx, UpredCI, 'b--','DisplayName','WNLM Fit & CI')
plot(0, Uis, 'b O', 'MarkerFaceColor', 'blue', 'DisplayName','Uis Endwall')
plot(0, Uis+Uis_bias, 'b O', 'DisplayName','Uis Bias')
plot(0, Uis-Uis_bias, 'b O', 'DisplayName','Uis Bias')
xlabel('Distance From Endwall [m]')
ylabel('Shock Velocity [m/s]')
hold off


ShockPassFrequency=Uis/WindowDiam; %[Hz] How fast the shock passes a test section window.
%Load and extract shock data
shock_file = strcat(file_base,num2str(date),'_',num2str(run_no,'%03d'),'_','ShockData.csv');
shock_data = readmatrix(shock_file,'NumHeaderLines',2);
time = shock_data(:,1);
PCB5 = shock_data(:,7);

%Load and extract dark/full data
dark_file = strcat(file_base,num2str(date),'_',num2str(run_no,'%03d'),'_','DarkFullData.csv');
dark_data = readmatrix(dark_file,'NumHeaderLines',2);

%Use File inputs for FROSH to calculate T5, P5, etc, and bias in T5 and P5 as a percent
EP_file = strcat(file_base,num2str(date),'_',num2str(run_no,'%03d'),'_','ExperimentParameters.csv');%Load and extract experimental parameters
EP_data = readmatrix(EP_file,'NumHeaderLines',2);%Load and extract experimental parameters
%Assign experimental parameter data to variables and assign bias values
T1 = EP_data(2,5)+273.15;
T1_bias = 5.0; %MDH change 11/9/2023.  RTD is 0.15K, but it's on the outside wall of the test section.
P1 = EP_data(1,5);
P1_bias = P1*0.0012;
solnMethod = 'EE';
nasaFile = 'THERM.DAT';
[Tout, Pout, rhoOut, gammaOut, VsoundOut, Mout, Urs, phi] = frosh('Uis',Uis,'T1',T1,'P1',P1,'testGasSpec',testGasSpec,'XtestGas',XtestGas,'solnMethod',solnMethod,'nasaFile',nasaFile);
T2=Tout(1); %[K]
P2=Pout(1); %[Pa]
T5=Tout(2); %[K]
P5=Pout(2); %[Pa]

timefromendwall=EndwalltoTestSection/Urs; %[s]
GasDynamicsFrequency=1/(timetoendwall+timefromendwall);
%This might need to go further down
UatIncident = U_fit.p1*EndwalltoTestSection + U_fit.p2; %[m/s] shock velocity at T_incident?
time_to_EW = roots([0.5*U_fit.p1 UatIncident -EndwalltoTestSection]); %[ms]
time_to_EW = time_to_EW(2)*1000; %convert s to ms

%% Pressure Trace
%Need to add code to predict and account for shock bifurcation.  Currently,
%buncha stuff starts in the middle of the bifurcation because Matlab thinks
%it's the ignition peak.

%Pressure Trace and IDT
pressurenoise=dark_data(:,2);
Kistler_Voltage = shock_data(:,2);%[V]
Pressure_atm=Kistler_Voltage.*(Kistler_psi_per_volt/14.7);% [atm]
Pressure_Pa=Kistler_Voltage.*Kistler_psi_per_volt*6894.76;% [pa]
Ti=Pressure_Pa.*T5/P5;%[K] This is a constant volume assumption. I don't think that's right. It's only close at state 5
disp('Temperature for CO2 absorbance calculations based on constant volume assumption. I **Think** it should be adiabatic, where PV^n (n=gamma) is constant, not volume..')
PXD_diameter=0.00508;%[m] = 0.2in - *approximate* diameter of kistler pxd
PXD_ShockPass_Frequency=Uis/PXD_diameter;
PXD_fft_Frequency=55000; %Prove to yourself this is a common noise frequency by uncommenting the next two lines of code.
%fftP=fft(pressurenoise);
%plot(OscopeSampleFrequency/length(pressurenoise)*(-length(pressurenoise)/2:length(pressurenoise)/2-1),abs(fftshift(fftP)),"LineWidth",3)
PXD_LP_Frequency=min([PXD_ShockPass_Frequency PXD_Amp_Frequency PXD_fft_Frequency]);
pressure_lp = lowpass(Kistler_Voltage,PXD_LP_Frequency,OscopeSampleFrequency); %y = lowpass(x,fpass,fs)
dPdt=gradient(pressure_lp,time);% point-by-point derivative of absorbCO2
PressureWindowCenter=int64(PretriggerIndex); %[index]center of "findpeaks" window.
ReflectedToIncidentBackwardsIndex=int64((timetoendwall+timefromendwall)*(OscopeSampleFrequency));
StartIncidentShockSearch=int64(PressureWindowCenter-(ReflectedToIncidentBackwardsIndex*2));
StopReflectedShockSearch=int64(PressureWindowCenter+3000);%Change 3000 to some a shock-specific number.
%Find incident and reflected shock peaks
[Ppks,Plocs,Pw,Pp] = findpeaks(dPdt(StartIncidentShockSearch:StopReflectedShockSearch),'MinPeakHeight',0.3*max(dPdt));%,'Npeaks',2,'MinPeakDistance',length(StartIncidentShockSearch:StopReflectedShockSearch)/2); %pks=peaks, locs=locations, w=width, p=prominence, 'MinPeakHeight' means ignore all peaks below 1500.
%[Ppks,Plocs,Pw,Pp] = findpeaks(dPdt(StartIncidentShockSearch:StopReflectedShockSearch),'MinPeakProminence',100,'Npeaks',2,'MinPeakDistance',((timetoendwall+timefromendwall)*OscopeSampleFrequency)/2); %pks=peaks, locs=locations, w=width, p=prominence, 'MinPeakHeight' means ignore all peaks below 1500.
Plocs=[Plocs(1)+StartIncidentShockSearch Plocs(2)+StartIncidentShockSearch];%redefine Plocs in terms of master index.
Ppks=[Ppks(1) Ppks(2)]; %make it same length as Plocs, in case there were more than 2 peaks.
PState1End=Plocs(1)-int64(((timetoendwall+timefromendwall)*OscopeSampleFrequency)/1);%index of State 1 start 
PState2Start=Plocs(1)+int64(((timetoendwall+timefromendwall)*OscopeSampleFrequency)/10); %index of State 2 start 
PState2End=Plocs(2)-int64(((timetoendwall+timefromendwall)*OscopeSampleFrequency)/3); %index of State 2 end
PState5start=Plocs(2)+int64(((timetoendwall+timefromendwall)*OscopeSampleFrequency)/5); %index of State 5 start
t_A_PXD = time(Plocs(2)); %but this won't work yet cuz your PXD filtering sucks.
%Find ignition from P
if strcmpi(UserInput{10},'y') %Reacting Tests
    [P_max,P_max_ind] = max(pressure_lp);
    [dPdtMax,dPdtMaxIndex] = max(dPdt(P_max_ind-3000:P_max_ind)); %Need to improve this too....point of maximum slope in the 3000 data points leading up to max P.
    t_ignPIndex=P_max_ind-3000+dPdtMaxIndex;
    t_ignP=time(t_ignPIndex);
    t_SW_P = time(Plocs(2));%[ms] %use max
    IDT_P=(t_ignP-t_SW_P)*1000; %[us]  
    IDT_P_bias='NA';  %not calculated YET.  working on it.
    PState5end=(t_ignPIndex); % state 5 ends when ignition begins. 
else %Non-reacting (Absorbing) tests
    Rarefaction_index = find(dPdt < 0);
    %State5end=Rarefaction_index(1);  %Rarefaction identification still not working right
    PState5end=PState5start+int64(((timetoendwall+timefromendwall)*OscopeSampleFrequency)*4);
    disp('State 5 start + State2 Lengthx4 used to determine Kistler State 5 end.  Still need to switch to using Rarefaction wave arrival')
    comments=strcat(comments,'Kistler-determined States unreliable with V1 code ');
   t_ignPIndex='NA';t_SW_P='NA'; IDT_P='NA';IDT_P_bias='NA';
end

figure(3)
hold on
yline(Pout(2)/101325,'DisplayName','FROSH P5')
plot(time,Pressure_atm,'DisplayName','PXD')
if strcmpi(UserInput{10},'y')
plot([t_ignP t_ignP], [0 max(Kistler_Voltage)*20/14.7],'b--','DisplayName','t Ign P')
else
end
plot([time(PState1End) time(PState1End)], [0 max(Kistler_Voltage)*20/14.7],'green--','DisplayName','State1End')
plot([time(PState2Start) time(PState2Start)], [0 max(Kistler_Voltage)*20/14.7],'k--','DisplayName','State2Start')
plot([time(PState2End) time(PState2End)], [0 max(Kistler_Voltage)*20/14.7],'k-','DisplayName','State2End')
plot([time(PState5start) time(PState5start)], [0 max(Kistler_Voltage)*20/14.7],'k:','DisplayName','State5start')
plot([time(PState5end) time(PState5end)], [0 max(Kistler_Voltage)*20/14.7],'k-.','DisplayName','State5End')
xlabel('time [ms]')
ylabel('Pressure [atm]')
legend
xlim([1.5 5])
ylim([0 max(Kistler_Voltage)*20/14.7])
hold off

figure(4)
hold on
plot(Kistler_Voltage,'DisplayName','1)pressure raw')
plot(pressure_lp,'DisplayName','2)pressure lp')
plot(dPdt./100,'DisplayName','4)dPdt./1000')
plot(Plocs,Ppks./100,'o','MarkerSize',12,'DisplayName','2)Inc & Ref Shock pks')
xlim([StartIncidentShockSearch StopReflectedShockSearch])
legend
hold off

%% HeNe
%Calculate absorbance time history and find reflected schlieren spike for HeNe laser. Can change which one is used by setting t_SW to the appropriate variable
if strcmpi(UserInput{9},'y')
%Shock File Data
HeNePitch = shock_data(:,4);
AvgHeNePitch=mean(HeNePitch);
StdHeNePitch=std(HeNePitch);
HeNeCatch = shock_data(:,5);
%Dark/Full File Data    
HeNePitchDark = dark_data(:,4);
AvgHeNePitchDark = mean(dark_data(:,4));
StdHeNePitchDark=std(dark_data(:,4));
HeNeCatchDark = dark_data(:,5);
AvgHeNeCatchDark = mean(dark_data(:,5));
StdHeNeCatchDark=std(dark_data(:,5));
absorbHeNe = -log((HeNeCatch./HeNePitch).*(HeNePitchDark./HeNeCatchDark));% With Common mode rejection
[~,spike_ind] = max(absorbHeNe(1:ceil(length(absorbHeNe)/2)));
t_A_HeNe = time(spike_ind); 
HeNedetectorFrequency=1000000; %[Hz] PDAVJ5 - DC-1MHz
absorbHeNe_lp = lowpass(absorbHeNe,HeNedetectorFrequency,OscopeSampleFrequency);
absorbHeNeBaseline=mean(absorbHeNe_lp(1:PretriggerIndex/2));
absorbHeNestd=std(absorbHeNe_lp(1:PretriggerIndex/2));
dHeNedt=gradient(absorbHeNe,time);% point-by-point derivative of absorbHeNe

if abs(-2*absorbHeNestd)<abs(absorbHeNeBaseline)||absorbHeNeBaseline>(2*absorbHeNestd) %If baseline absorbance is far from zero, correct it...
   absorbHeNe_lp=absorbHeNe-absorbHeNeBaseline; %Correct to zero so you can still use this to find State 2 & 5 peaks, but:
else  %If baseline absorbance is far from zero, don't use absorbance data, but use the schlieren spikes to ID states 2 and 5.  
end
%dHeNedtcontainsImaginary=~isreal(dHeNedt);
%if dHeNedtcontainsImaginary
%    disp('dHeNe/dt contains imaginary numbers. The signal was probably bad, and all your HeNe-derived data will be written as NA');
%    State1_HeNe_abs='NA'; State2_HeNe_abs='NA'; Std_State2_HeNe_abs='NA'; State5_HeNe_abs='NA'; Std_State5_HeNe_abs='NA'; HeNePitchVariation='NA';
%else
[HeNepks,HeNelocs,HeNew,HeNep] = findpeaks(absorbHeNe_lp(1:StopReflectedShockSearch),'MinPeakHeight',2*max(absorbHeNe_lp(100:PretriggerIndex/1.5)),'MinPeakDistance',((timetoendwall+timefromendwall)*OscopeSampleFrequency)/2); %pks=peaks, locs=locations, w=width, p=prominence, 'MinPeakHeight' means ignore all peaks below 1500.
HeNeState1End=HeNelocs(1)-int64(((timetoendwall+timefromendwall)*OscopeSampleFrequency)/1);%index of State 1 start 
HeNeState2Start=HeNelocs(1)+int64(((timetoendwall+timefromendwall)*OscopeSampleFrequency)/10); %index of State 2 start 
HeNeState2End=HeNelocs(find(HeNepks==max(HeNepks)))-int64(((timetoendwall+timefromendwall)*OscopeSampleFrequency)/3); %index of State 2 end
HeNeCatch1= mean(shock_data(3:HeNeState1End,5));
HeNeCatchVariation=HeNeCatchDark/HeNeCatch1;
HeNePitchVariation=AvgHeNePitchDark/AvgHeNePitch;
if strcmpi(UserInput{10},'y') %Reacting Tests
    HeNeState5end=(t_ignPIndex); % state 5 ends when ignition begins. 
    disp('PXD-based Ignition time used to determine HeNe State 5 end')
    comments=strcat(comments,'PXD-based tign used to determine HeNe State 5 end');
else %Non-reacting (Absorbing) tests
    Rarefaction_index = find(dPdt < 0);
    %State5end=Rarefaction_index(1);  %Rarefaction identification still not working right
    HeNeState5end=HeNeState5start+int64(((timetoendwall+timefromendwall)*OscopeSampleFrequency)*4);
    disp('State 5 start + 0.5ms used to determine HeNe State 5 end.  Still need to switch to using Rarefaction wave arrival')
    comments=strcat(comments,'State 5 start + 0.5ms used to determine HeNe State 5 end');
end

HeNeState5start=HeNelocs(find(HeNepks==max(HeNepks)))+int64(((timetoendwall+timefromendwall)*OscopeSampleFrequency)/5); %index of State 5 start

if abs(-2*absorbHeNestd)<abs(absorbHeNeBaseline)||absorbHeNeBaseline>(2*absorbHeNestd) %If baseline absorbance is far from zero, don't use the absorbance data.
    State1_HeNe_abs='NA'; State2_HeNe_abs='NA'; Std_State2_HeNe_abs='NA'; State5_HeNe_abs='NA'; Std_State5_HeNe_abs='NA';%Don't use any of the absorbance data
    State2_sigma_HC='NA'; State2_sigma_HC_bias='NA';  State5_sigma_HC='NA'; State5_sigma_HC_bias='NA'; X_HC ='NA'; X_HC_bias='NA';  
    disp('Baseline HeNe absorbance 2*std away from 0. HeNeAbs NA.')
    comments=strcat(comments,'BL Abs_HeNe > 2*std away from 0. ');
else
 State1_HeNe_abs=mean(absorbHeNe(1:HeNeState1End));
    State2_HeNe_abs=mean(absorbHeNe(HeNeState2Start:HeNeState2End));
    Std_State2_HeNe_abs=std(absorbHeNe(HeNeState2Start:HeNeState2End));
    State5_HeNE_abs=mean(absorbHeNe(HeNeState5start:HeNeState5end));
    Std_State5_HeNe_abs=std(absorbHeNe(HeNeState5start:HeNeState5end));
    State2_sigma_HC=(State2_HeNe_abs*Ru*T2)/(XtestGas(1)*P2*TubeDiameter);  %not calculated YET.  working on it.
    State2_sigma_HC_bias='NA';  %not calculated YET.  working on it.
    State5_sigma_HC=(State2_HeNe_abs*Ru*T5)/(XtestGas(1)*P5*TubeDiameter);  
    State5_sigma_HC_bias='NA';  %not calculated YET.  working on it.
    X_HC ='NA';  %not calculated YET.  working on it.
    X_HC_bias='NA';  %not calculated YET.  working on it..
end


figure(5)
hold on
plot(time, HeNePitch./mean(HeNePitch(1:PretriggerIndex/2)),'r-','DisplayName','HeNe Ref Norm')
plot(time, HeNeCatch./HeNeCatch1,'b-','DisplayName','HeNe It Norm')
%plot(time,HeNeCatchDark./HeNeCatch1,'b-.','DisplayName','HeNeI0 Norm') %This one in particular makes Matlab choke.
plot(time, absorbHeNe,'k-','DisplayName','HeNe absorbance')
plot([time(HeNeState1End) time(HeNeState1End)], [0 1],'green--','DisplayName','State1End')
plot([time(HeNeState2Start) time(HeNeState2Start)], [0 1],'k--','DisplayName','State2Start')
plot([time(HeNeState2End) time(HeNeState2End)], [0 1],'k-','DisplayName','State2End')
plot([time(HeNeState5start) time(HeNeState5start)], [0 1],'k:','DisplayName','State5start')
plot([time(HeNeState5end) time(HeNeState5end)], [0 1],'k-.','DisplayName','State5End')
xlabel('time [ms]')
ylabel('Signal [-]')
legend
xlim([1.5 5])
hold off

figure(6)
%plot(dHeNedt./1000,'DisplayName','4)dHeNedt./1000')
plot(absorbHeNe,'DisplayName','1)absorbHeNe raw')
hold on
plot(absorbHeNe_lp,'DisplayName','2)absorbHeNe lp')
yline(absorbHeNeBaseline,'DisplayName',['Real Baseline Abs'])
plot(HeNelocs,HeNepks,'o','MarkerSize',12,'DisplayName','2)Inc & Ref Shock pks')
legend
hold off
t_SW = t_A_HeNe; %[ms] time of REFLECTED? shock arrival at sidewall.
T_incident = time(HeNelocs(1));%[ms] time of INCIDENT shock arrival at ENDWALL?
%Uat incident code goes here if it doesn't work above.
t_EW = T_incident + time_to_EW;%[ms] time of incident shock arrival at end wall
%end %end of "if contains imaginary' statement
else
    State2_HeNe_abs='NA';Std_State2_HeNe_abs='NA';State5_HeNe_Abs='NA';Std_State2_HeNe_abs='NA';AvgHeNePitchDark='NA';StdHeNePitchDark='NA';
    AvgHeNeCatchDark='NA';StdHeNeCatchDark='NA'; AvgHeNePitch='NA'; StdHeNePitch='NA'; State1_HeNe_abs='NA'; Std_State5_HeNe_abs='NA';HeNePitchVariation='NA'; ...%HeNe Laser Outputs
    State2_sigma_HC='NA'; State2_sigma_HC_bias='NA'; X_HC='NA'; X_HC_bias='NA'; 
    disp('No HeNe Data used. t_SW comes from....?' )
    t_SW = 'NA figure out t_A_PXD';
    T_incident = 'NA figure out t_A_PXD';%[ms] time of INCIDENT shock arrival at ENDWALL?
    t_EW ='NA figure out t_A_PXD';
end


%%Ask Michael about details!
%These are all calculated using the HeNe code block, but the timings COULD
%also be derived from the P and CO2 traces.  Determine uncertainty associated
%with each of the three methods, choose method with least uncertainty, and
%report it.


%Calculate time of incident shock arrival at PCB5
UatPCB5 = U_fit.p1*EndwalltoPCB5 + U_fit.p2;
PCB5_lp = lowpass(PCB5,1000,10e6);
t_PCB5 = time(find(PCB5_lp>0.002,1));

%Frosh Bias [% % % % K atm]
[T2_bias P2_bias T5_bias, P5_bias, T5, P5_atm] = frosh_bias_V2('Uis',Uis,'Uis_bias',Uis_bias,'T1',T1,'T1_bias',T1_bias,'P1',P1,'P1_bias',P1_bias,'testGasSpec',testGasSpec,'XtestGas',XtestGas,'XtestGas_bias',XtestGas_bias,'solnMethod',solnMethod,'nasaFile',nasaFile);

%OH* Side Wall IDT
if strcmpi(UserInput{7},'y') && strcmpi(UserInput{9},'y')%if OH Side Wall data {6} was taken &&(Still need HeNe data {9} specifically to determine IDT.)
%Check for signal noise frequencies by uncommenting the following three lines.
%OH_SW_noise = dark_data(:,3);
%fftOHSW=fft(OH_SW_noise);
%plot(OscopeSampleFrequency/length(OH_SW)*(-length(OH_SW)/2:length(OH_SW)/2-1),abs(fftshift(fftOHSW)),"LineWidth",3)
OH_SW = shock_data(:,3); 
OH_SW_max=max(OH_SW);
OH_SW_lp = lowpass(OH_SW,PhotodetectorFrequency,OscopeSampleFrequency);
[OH_SW_lp_max,OH_SW_max_ind] = max(OH_SW_lp);
OH_SW_sigmoid=fit((PretriggerIndex:OH_SW_max_ind)',OH_SW_lp(PretriggerIndex:OH_SW_max_ind),'logistic');
OH_SW_sigmoid_coefficients=coeffvalues(OH_SW_sigmoid);
if OH_SW_sigmoid_coefficients(1)>2*OH_SW_lp_max||OH_SW_sigmoid_coefficients(1)<0.5*OH_SW_lp_max %If the sigmoid fit is really bad, just use old way of finding max slope.
dEdt_OH_SW=gradient(OH_SW_lp,time); 
disp('OH SW sigmoid fit sucked.  dEdt comes from gradient(OH_SW_lp)')
comments=strcat(comments,'IDT OH SW from Gradient(OH_SW_lp)');
else
dEdt_OH_SW=gradient(OH_SW_sigmoid(1:length(OH_SW)),time); %otherwise use Sigmoid slope.
end
PreshockE_OH_SW=mean(OH_SW_lp(1:spike_ind)); %Pre-ignition baseline
[dEdt_OH_SWpks,dEdt_OH_SWlocs,dEdt_OH_SWw,dEdt_OH_SWp] = findpeaks(dEdt_OH_SW,'MinPeakHeight',0.9*max(dEdt_OH_SW)); %pks=peaks, locs=locations, w=width, p=prominence, 'MinPeakHeight' means ignore all peaks below 1500.
EmissionIntercept_OH_SW=OH_SW_sigmoid(dEdt_OH_SWlocs(1))-(dEdt_OH_SWpks(1)*time(dEdt_OH_SWlocs(1))); %Y-intercept for line through dEdtMAx
dEdtMaxLine_OH_SW=dEdt_OH_SWpks(1).*time+EmissionIntercept_OH_SW; %Full equation for line through dEdtMAx
[~, index_intercept_OH_SW] = min(abs(dEdtMaxLine_OH_SW-PreshockE_OH_SW));%Det
t_ignOH_SW=time(index_intercept_OH_SW);
IDT_OH_SW=(t_ignOH_SW-t_SW)*1000; %[us]
IDT_OH_SW_bias='NA';  %not calculated YET.  working on it.

figure(7) 
hold on
plot(time,OH_SW,'DisplayName','OH SW')
plot(time,OH_SW_lp,'DisplayName','OH SW lp')
plot(time,OH_SW_sigmoid(1:length(OH_SW)),'DisplayName','OH sigmoid fit')
%plot([T_incident T_incident],[0 max(OH_SW)],'k','DisplayName','t incident')
plot([t_SW t_SW],[0 max(OH_SW)],'r--','DisplayName','t SW')
plot([time(dEdt_OH_SWlocs(1)) time(dEdt_OH_SWlocs(1))], [0 max(OH_SW)],'k-.','DisplayName','dEdt Max')
plot(time,dEdtMaxLine_OH_SW)
yline(PreshockE_OH_SW,'DisplayName','Baseline Emission')
plot([t_ignOH_SW t_ignOH_SW], [0 max(OH_SW)],'r-.','DisplayName','t Ign SW')
xlim([1.5 5])
ylim([0 1.5*max(OH_SW)])
xlabel('time [ms]')
ylabel('Volts [V]')
legend
hold off

else
    IDT_OH_SW='NA';t_ignOH_SW='NA';IDT_OH_SW_bias='NA';
    disp('No OH SW Data used')
end


%OH* End Wall IDT
if strcmpi(UserInput{6},'y') && strcmpi(UserInput{9},'y')%if OH End Wall data {6} was taken &&(Still need HeNe data {9} specifically to determine IDT.)
OH_EW = -shock_data(:,6);
OH_EW_lp = lowpass(OH_EW,1e5,OscopeSampleFrequency);
[OH_EW_lp_max,OH_EW_max_ind] = max(OH_EW_lp);
OH_EW_sigmoid=fit((PretriggerIndex:OH_EW_max_ind)',OH_EW_lp(PretriggerIndex:OH_EW_max_ind),'logistic');
OH_EW_sigmoid_coefficients=coeffvalues(OH_EW_sigmoid);
if OH_EW_sigmoid_coefficients(1)>2*OH_EW_lp_max||OH_EW_sigmoid_coefficients(1)<0.5*OH_EW_lp_max %If the sigmoid fit is really bad, just use old way of finding max slope.
dEdt_OH_EW=gradient(OH_EW_lp,time); %
disp('OH EW sigmoid fit sucked.  dEdt comes from gradient(OH_EW_lp)')
comments=strcat(comments,'IDT OH EW from Gradient(OH_EW_lp)');
else
dEdt_OH_EW=gradient(OH_EW_sigmoid(1:length(OH_EW)),time); %otherwise use Sigmoid slope.
end
PreshockE_OH_EW=mean(OH_EW_lp(1:spike_ind)); %Pre-ignition baseline
[dEdt_OH_EWpks,dEdt_OH_EWlocs,dEdt_OH_EWw,dEdt_OH_EWp] = findpeaks(dEdt_OH_EW,'MinPeakHeight',0.1*max(dEdt_OH_EW)); %pks=peaks, locs=locations, w=width, p=prominence, 'MinPeakHeight' means ignore all peaks below 1500.
EmissionIntercept_OH_EW=OH_EW_sigmoid(dEdt_OH_EWlocs(1))-(dEdt_OH_EWpks(1)*time(dEdt_OH_EWlocs(1))); %Y-intercept for line through dEdtMAx
dEdtMaxLine_OH_EW=dEdt_OH_EWpks(1).*time+EmissionIntercept_OH_EW; %Full equation for line through dEdtMAx
[~, index_intercept_OH_EW] = min(abs(dEdtMaxLine_OH_EW-PreshockE_OH_EW));%Det
t_ignOH_EW=time(index_intercept_OH_EW);
IDT_OH_EW=(t_ignOH_EW-t_EW)*1000; %[us]
IDT_OH_EW_bias='NA';  %not calculated YET.  working on it.

figure(8)
hold on
plot(time,(OH_EW),'DisplayName','OH EW')
plot(time,(OH_EW_lp),'DisplayName','OH EW lp')
plot(time,OH_EW_sigmoid(1:length(OH_EW)),'DisplayName','OH sigmoid fit')
plot([t_EW t_EW],[0 max(OH_EW)],'b--','DisplayName','t EW')
plot([t_ignOH_EW t_ignOH_EW], [0 max(OH_EW)],'b-.','DisplayName','t Ign EW')
plot([time(dEdt_OH_EWlocs(1)) time(dEdt_OH_EWlocs(1))], [0 max(OH_EW)],'k-.','DisplayName','dEdt Max')
plot(time,dEdtMaxLine_OH_EW)
yline(PreshockE_OH_EW,'DisplayName','Baseline Emission')
xlim([1.5 5])
ylim([0 1.5*max(OH_EW)])
xlabel('time [ms]')
ylabel('Volts [V]')
legend
hold off

else
    IDT_OH_EW='NA';t_ignOH_EW='NA';IDT_OH_EW_bias='NA';
    disp('No OH EW Data used')
end


%CO2 Laser Absorption Calculations
if strcmpi(UserInput{8},'y')
%Shock Data File
CO2Pitch = shock_data(:,8); %Iref-signal (dark offset applied)
AvgCO2Pitch=mean(CO2Pitch);
StdCO2Pitch=std(CO2Pitch);
CO2Catch = shock_data(:,9);%It-signal (dark offset applied)
%Dark/Full File Data
CO2PitchDark = mean(dark_data(:,8)); %Iref-vaccuum (dark offset applied)
StdCO2PitchDark=std(dark_data(:,8)); 
CO2CatchDark = mean(dark_data(:,9));%I0-vacuum (dark offset applied)
StdCO2CatchDark=std(dark_data(:,9));
I0CO2_bias=StdCO2CatchDark; %Claims bias is just std dev of signal
It_CO2=CO2Catch.*(CO2PitchDark./CO2CatchDark); %It CMR-corrected
absorbCO2 = -log((CO2Catch./CO2Pitch).*(CO2PitchDark./CO2CatchDark));% With Common mode rejection (Dark signal offset applied in LabView)
CO2detectorFrequency=1000000; %Hz
absorbCO2_lp = lowpass(absorbCO2,GasDynamicsFrequency,OscopeSampleFrequency);
dCO2dt=gradient(absorbCO2_lp,time);% point-by-point derivative of absorbCO2
dCO2dtcontainsImaginary=~isreal(dCO2dt);
absorbCO2Baseline=mean(absorbCO2_lp(1:PretriggerIndex/2));
absorbCO2std=mean(absorbCO2_lp(1:PretriggerIndex/2));

if abs(-2*absorbCO2std)<abs(absorbCO2Baseline)||absorbCO2Baseline>(2*absorbCO2std) %If baseline absorbance is far from zero, don't use the absorbance data.
    absorbCO2=absorbCO2_lp-absorbCO2Baseline;
    disp('Baseline CO2 absorbance was < 0. Signal might be unreliable')
    comments=strcat(comments,'Baseline CO2 absorbance was < 0. ');
else %Do nothing.
end

%if dCO2dtcontainsImaginary
%    disp('dCO2/dt contains imaginary numbers. The signal was probably bad, and all your CO2-derived data will be written as NA');
%    State1_CO2_abs='NA'; State2_CO2_abs='NA'; Std_State2_CO2_abs='NA'; State5_CO2_abs='NA'; Std_State5_CO2_abs='NA'; CO2PitchVariation='NA';
%else
%[CO2pks,CO2locs,CO2w,CO2p] = findpeaks(absorbCO2_lp,'MinPeakHeight',3*max(absorbCO2_lp(100:PretriggerIndex/2)),'MinPeakDistance',((timetoendwall+timefromendwall)*OscopeSampleFrequency)/2); %pks=peaks, locs=locations, w=width, p=prominence, 'MinPeakHeight' means ignore all peaks below 1500.
[CO2pks,CO2locs,CO2w,CO2p] = findpeaks(dCO2dt,'MinPeakHeight',3*max(dCO2dt(100:PretriggerIndex/2)),'MinPeakDistance',((timetoendwall+timefromendwall)*OscopeSampleFrequency)/2); %pks=peaks, locs=locations, w=width, p=prominence, 'MinPeakHeight' means ignore all peaks below 1500.
%[CO2pks,CO2locs,CO2w,CO2p] = findpeaks(dCO2dt,'MinPeakHeight',0.25*max(dCO2dt),'MinPeakDistance',((timetoendwall+timefromendwall)*OscopeSampleFrequency)/2,'MinPeakWidth',0.5*max(CO2w)); %pks=peaks, locs=locations, w=width, p=prominence, 'MinPeakHeight' means ignore all peaks below 1500.
State1End=CO2locs(1)-int64(((timetoendwall+timefromendwall)*OscopeSampleFrequency)/1);%index of State 1 start 
State2Start=CO2locs(1)+int64(((timetoendwall+timefromendwall)*OscopeSampleFrequency)/10); %index of State 2 start 
State2End=CO2locs(find(CO2pks==max(CO2pks)))-int64(((timetoendwall+timefromendwall)*OscopeSampleFrequency)/3); %index of State 2 end based on Reflected shock Schlieren
CO2Catch1= mean(shock_data(3:State1End,9));
CO2CatchVariation=CO2CatchDark/CO2Catch1;
CO2PitchVariation=CO2PitchDark/AvgCO2Pitch;
State5start=CO2locs(find(CO2pks==max(CO2pks)))+int64(((timetoendwall+timefromendwall)*OscopeSampleFrequency)/5); %index of State 5 start based on Reflected shock Schlieren
if strcmpi(UserInput{10},'y') %Reacting Tests:  state 5 ends when ignition begins.
    if strcmpi(UserInput{6},'y')
         State5end=(index_intercept_OH_EW)-500; % If OH endwall was used, use that, cuz it should happen earliest    
         disp('OH Endwall tign used for CO2 State 5 end')
         comments=strcat(comments,'OH EW tign used for CO2 State 5 end');
    else 
        if strcmpi(UserInput{7},'y') %if OH endwall wasn't used, but OH sidewall was, 
           State5end=(index_intercept_OH_SW)-500; % Use OH sidewall, cuz it's more reliable than pressure
           disp('OH Sidewall tign used for CO2 State 5 end')
           comments=strcat(comments,'OH SW tign used for CO2 State 5 end');
        else
            State5end=(t_ignPIndex)-500; % If you have not other choice, use P
            disp('PXD tign used for CO2 State 5 end')
            comments=strcat(comments,'PXD tign used for CO2 State 5 end');
        end
    end
else %Non-reacting (Absorbing) tests
    Rarefaction_index = find(dPdt < 0);
    %State5end=Rarefaction_index(1);  %Rarefaction identification still not working right
    State5end=State5start+int64(((timetoendwall+timefromendwall)*OscopeSampleFrequency)*4);
    disp('State 5 start State2 Lengthx4  used to determine CO2 Laser State 5 end.  Still need to switch to using Rarefaction wave arrival')
    comments=strcat(comments,'State 5 start State2 Lengthx4 used to determine CO2 State 5 end');
    
end

if abs(-2*absorbCO2std)<abs(absorbCO2Baseline)||absorbCO2Baseline>(2*absorbCO2std) %If baseline absorbance is far from zero, don't use the absorbance data.
State1_CO2_abs='NA'; State2_CO2_abs='NA';Std_State2_CO2_abs='NA';State5_CO2_abs='NA';Std_State5_CO2_abs='NA';State2sigma_CO2='NA';State2sigma_CO2_bias='NA'; State2_X_CO2='NA';  State2_X_CO2_bias ='NA';  
else
State1_CO2_abs=mean(absorbCO2(1:State1End));
State2_CO2_abs=mean(absorbCO2(State2Start:State2End));
State2_It_CO2=mean(It_CO2(State2Start:State2End));
State2_It_CO2_bias=std(It_CO2(State2Start:State2End));%Claims bias is just std dev of signal
Std_State2_CO2_abs=std(absorbCO2(State2Start:State2End)); 
State5_CO2_abs=mean(absorbCO2(State5start:State5end));
Std_State5_CO2_abs=std(absorbCO2(State5start:State5end));
State5_It_CO2=mean(It_CO2(State5start:State5end));
State5_It_CO2_bias=std(It_CO2(State5start:State5end));%Claims bias is just std dev of signal
if strcmpi(testGasSpec{1},'C2H4')
State2sigma_CO2=(State2_CO2_abs*Ru*T2)/(XtestGas(1)*P2*TubeDiameter);  
State2_delsigma_delx=(-State2_CO2_abs*Ru*T2)/(XtestGas(1)^2*TubeDiameter*P2); %influence coefficient - X_C2H4
State2_delsigma_delP=(-State2_CO2_abs*Ru*T2)/(XtestGas(1)*TubeDiameter*P2^2); %influence coefficient - P2 
State2_delsigma_delL=(-State2_CO2_abs*Ru*T2)/(XtestGas(1)*TubeDiameter^2*P2); %influence coefficient - L
State2_delsigma_delT=(-State2_CO2_abs*Ru)/(XtestGas(1)*TubeDiameter*P2); %influence coefficient - T2
State2_delsigma_delIt=(-Ru*T2)/(XtestGas(1)*TubeDiameter*P2*State2_It_CO2); %influence coefficient - It
State2_delsigma_delI0=(-Ru*T2)/(XtestGas(1)*TubeDiameter*P2*CO2CatchDark); %influence coefficient - I0
State2sigma_CO2_bias=sqrt((State2_delsigma_delx*XtestGas_bias(1))^2+(State2_delsigma_delP*P2_bias)^2+(State2_delsigma_delL*Bias_TubeDiameter)^2+(State2_delsigma_delT*T2_bias)^2+(State2_delsigma_delIt*State2_It_CO2_bias)^2+(State2_delsigma_delI0*I0CO2_bias)^2);  %[m^2/mol]
State5sigma_CO2=(State2_CO2_abs*Ru*T5)/(XtestGas(1)*P5*TubeDiameter);
State5_delsigma_delx=(-State5_CO2_abs*Ru*T5)/(XtestGas(1)^2*TubeDiameter*P5); %influence coefficient - X_C2H4
State5_delsigma_delP=(-State5_CO2_abs*Ru*T5)/(XtestGas(1)*TubeDiameter*P5^2); %influence coefficient - P2 
State5_delsigma_delL=(-State5_CO2_abs*Ru*T5)/(XtestGas(1)*TubeDiameter^2*P5); %influence coefficient - L
State5_delsigma_delT=(-State5_CO2_abs*Ru)/(XtestGas(1)*TubeDiameter*P5); %influence coefficient - T2
State5_delsigma_delIt=(-Ru*T5)/(XtestGas(1)*TubeDiameter*P5*State5_It_CO2); %influence coefficient - It
State5_delsigma_delI0=(-Ru*T5)/(XtestGas(1)*TubeDiameter*P5*CO2CatchDark); %influence coefficient - I0
State5sigma_CO2_bias=sqrt((State5_delsigma_delx*XtestGas_bias(1))^2+(State5_delsigma_delP*P5_bias)^2+(State5_delsigma_delL*Bias_TubeDiameter)^2+(State5_delsigma_delT*T5_bias)^2+(State2_delsigma_delIt*State5_It_CO2_bias)^2+(State5_delsigma_delI0*I0CO2_bias)^2);  %[m^2/mol]
else
State2sigma_CO2='NA';State2sigma_CO2_bias='NA';State5sigma_CO2='NA'; State5sigma_CO2_bias='NA'; 
end
State2_Sigma_Pinkowski=5.97+44.04*exp(-((T2-370.56)/259.77))+42.3*exp(-((T2-370.56)/259.792)); %copied from my PC. Double-check all coefficients with Pink. paper
State5_Sigma_Pinkowski=5.97+44.04*exp(-((T5-370.56)/259.77))+42.3*exp(-((T5-370.56)/259.792)); %copied from my PC. Double-check all coefficients with Pink. paper
Statei_Sigma_Pinkowski=5.97+44.04*exp(-((Ti-370.56)/259.77))+42.3*exp(-((Ti-370.56)/259.792)); %copied from my PC. Double-check all coefficients with Pink. paper
State2_X_CO2=(State2_CO2_abs*Ru*T2)/(State2_Sigma_Pinkowski*P2*TubeDiameter);  %Single point, average value
State2_X_Error=1-State2_X_CO2/XtestGas(1); %percent error compared to manometric measurement
State2_X_CO2_bias ='NA';  %not calculated YET.  working on it.
State5_X_CO2=(State5_CO2_abs*Ru*T5)/(State5_Sigma_Pinkowski*P5*TubeDiameter);
State5_X_Error=1-State5_X_CO2/XtestGas(1); %percent error compared to manometric measurement
State5_X_CO2_bias ='NA';  %not calculated YET.  working on it.
Statei_X_CO2=(absorbCO2_lp.*Ru.*Ti)./(Statei_Sigma_Pinkowski.*Pressure_Pa.*TubeDiameter);  %calculatated X at every point

figure(09)
plot(Statei_X_CO2,'DisplayName','X CO_2 raw')
hold on
yline(State2_X_CO2,'green--','DisplayName','State2 X CO_2')
yline(State5_X_CO2,'r-','DisplayName','State5 X CO_2')
yline(XtestGas(1),'DisplayName','Manometric X CO_2')
legend
hold off
end


figure(10)
hold on
plot(time, CO2Pitch./mean(CO2Pitch(1:15000)),'r-','DisplayName','CO2 Ref Norm')
plot(time, CO2Catch./CO2Catch1,'b-','DisplayName','CO2 It Norm')
%plot(time,CO2CatchDark./CO2Catch1,'b-.','DisplayName','CO2I0 Norm') %This one in particular makes Matlab choke.
plot(time, absorbCO2,'k-','DisplayName','CO2 absorbance')
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

figure(11)
plot(absorbCO2,'DisplayName','1)absorbCO_2 raw')
hold on
plot(dCO2dt./100,'DisplayName','4)dCO_2dt./1000')
plot(CO2locs,CO2pks./100,'o','MarkerSize',12,'DisplayName','2)Inc & Ref Shock pks')
legend
hold off
%end %end of imaginary number if statement.

else 
    disp('No CO2 Laser Data used')
    State1_CO2_abs='NA'; State2_CO2_abs='NA'; Std_State2_CO2_abs='NA'; State5_CO2_abs='NA'; Std_State5_CO2_abs='NA'; CO2PitchVariation='NA'; CO2PitchDark='NA';
    StdCO2PitchDark='NA'; CO2CatchDark='NA'; StdCO2CatchDark='NA'; AvgCO2Pitch='NA'; StdCO2Pitch='NA'; State2sigma_CO2='NA'; State2sigma_CO2_bias='NA'; State2_X_CO2='NA'; State2_X_CO2_bias='NA'; CO2CatchVariation='NA';
end 
%Wait to bring up comments dialog until after Figure 5 is closed.
disp('NO DATA WRITTEN YET! Script is paused so you can look at graphs.  Analyze the graphs, then press ENTER and add any comments.')
pause

comments = inputdlg({'Any comments to add?  Dont use commas.'},'User Comments',[5 50],{comments});
%Single-point data outputs.  Written to a text file in the data processing folder
T=[date {num2str(run_no,'%03d')} testGasSpec XtestGas phi T1 P1 Uis Uis_bias Urs ... %Initital Conditions and Frosh Input
    P2 P2_bias T2 T2_bias gammaOut(2) P5 P5_bias T5 T5_bias gammaOut(3) 1000/T5 ... %Frosh Output.  States 2 and 5
    IDT_OH_EW IDT_OH_EW_bias t_EW t_ignOH_EW IDT_OH_SW IDT_OH_SW_bias t_SW t_ignOH_SW IDT_P IDT_P_bias ...%Ignition Delay Time from all 3 sources
    CO2PitchDark StdCO2PitchDark CO2CatchDark StdCO2CatchDark AvgCO2Pitch StdCO2Pitch State1_CO2_abs State2_CO2_abs Std_State2_CO2_abs State5_CO2_abs Std_State5_CO2_abs CO2PitchVariation...%CO2 laser Outputs
    State2sigma_CO2 State2sigma_CO2_bias State2_X_CO2 State2_X_CO2_bias CO2CatchVariation...%C2H4 data calculated from CO2
    AvgHeNePitchDark StdHeNePitchDark AvgHeNeCatchDark StdHeNeCatchDark AvgHeNePitch StdHeNePitch State1_HeNe_abs State2_HeNe_abs Std_State2_HeNe_abs Std_State5_HeNe_abs HeNePitchVariation ...%HeNe Laser Outputs
    State2_sigma_HC State2_sigma_HC_bias X_HC X_HC_bias ...%HC data calculated from HeNe
    comments]; %any user comments to alert the user to bad data or explanations of unusual data.
writecell(T,'CST_ProcessorOutput.txt','Delimiter',',','WriteMode','Append')
disp('Data written. All done.')