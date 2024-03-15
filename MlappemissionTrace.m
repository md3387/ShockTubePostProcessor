function [emission, emissionLP, t_ign, IDT, IDTBias,comments]=MlappemissionTrace(name,t_arrival, Pretrigger, time,emission ,LowPassFrequency,OscopeSampleFrequency, spike_ind,comments)

%Check for signal noise frequencies by uncommenting the following three lines.
%OH_SW_noise = dark_data(:,3);
%fftOHSW=fft(OH_SW_noise);
%plot(OscopeSampleFrequency/length(OH_SW)*(-length(OH_SW)/2:length(OH_SW)/2-1),abs(fftshift(fftOHSW)),"LineWidth",3)
%Pretrigger=str2double(UserInput{11}); %[ms]
PretriggerIndex=(Pretrigger/1000)*(OscopeSampleFrequency);
%time = shock_columns(:,1);
%emission = shock_columns(:,2);
emissionSG= sgolayfilt(emission,2,51);
BWorder = 2; % Order of the Butterworth filter
[b, a] = butter(BWorder, LowPassFrequency/ (OscopeSampleFrequency/ 2), 'low');% Design Butterworth filter
emissionBW = filtfilt(b, a, emissionSG);% Apply filter using filtfilt to achieve zero-phase filtering
emissionMerge = sqrt(abs(emissionBW .* emissionSG));

emissionLP = emissionMerge;%lowpass(emission,LowPassFrequency,OscopeSampleFrequency);
%[emissionLPmax,emissionLPmaxInd] = max(emissionLP);
%emissionSigmoid=fit((PretriggerIndex:emissionLPmaxInd)',emissionLP(PretriggerIndex:emissionLPmaxInd),'logistic');
%emissionSigmoidCoefficients=coeffvalues(emissionSigmoid);
%if emissionSigmoidCoefficients(1)>2*emissionLPmax||emissionSigmoidCoefficients(1)<0.5*emissionLPmax %If the sigmoid fit is really bad, just use old way of finding max slope.
dEdt=gradient(emissionLP,time); %
%disp(strcat(name,'sigmoid fit sucked.  dEdt comes from LowPass gradient'))
%comments=strcat(comments,' ',name,' IDT from LowPass gradient');
%else
%dEdt=gradient(emissionSigmoid(1:length(emission)),time); %otherwise use Sigmoid slope.
%end
PreshockEmission=mean(emissionLP(1:spike_ind)); %Pre-ignition baseline
[dEdtpks,dEdtlocs,dEdtw,dEdtp] = findpeaks(dEdt,'MinPeakHeight',0.9*max(dEdt)); %pks=peaks, locs=locations, w=width, p=prominence, 'MinPeakHeight' means ignore all peaks below 1500.
%EmissionIntercept=emissionSigmoid(dEdtlocs(1))-(dEdtpks(1)*time(dEdtlocs(1))); %Y-intercept for line through dEdtMAx
EmissionIntercept=emissionLP(dEdtlocs(1))-(dEdtpks(1)*time(dEdtlocs(1))); %Y-intercept for line through dEdtMAx
dEdtMaxLine=dEdtpks(1).*time+EmissionIntercept; %Full equation for line through dEdtMAx
[~, index_intercept] = min(abs(dEdtMaxLine-PreshockEmission));%Det
t_ign=time(index_intercept);
IDT=(t_ign-t_arrival)*1000; %[us]
IDTBias='NA';  %not calculated YET.  working on it.

%{
figureIndex=figureIndex+1;
figure(figureIndex)
hold on
plot(time,emission,'DisplayName',name)
plot(time,emissionLP,'DisplayName',strcat(name,'LP'))
%plot(time,emissionSigmoid(1:length(emission)),'DisplayName','sigmoid fit')
plot([t_arrival t_arrival],[0 max(emission)],'r--','DisplayName','t_arrival')
plot([time(dEdtlocs(1)) time(dEdtlocs(1))], [0 max(emission)],'k-.','DisplayName','dEdt Max')
plot(time,dEdtMaxLine)
yline(PreshockEmission,'DisplayName','Baseline Emission')
plot([t_ign t_ign], [0 max(emission)],'r-.','DisplayName','t Ign')
xlim([1.5 5])
ylim([0 1.5*max(emission)])
xlabel('time [ms]')
ylabel('Volts [V]')
legend
hold off
%}