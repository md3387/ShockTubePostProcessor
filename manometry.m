function [XtestGas, XtestGas_bias, MW] = manometry(Pmix,testGasSpec)

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