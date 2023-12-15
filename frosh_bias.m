function [B_T2 B_P2 B_T5, B_P5, T5, P5] = frosh_bias_V2(varargin)

global XtestGas_temp Num_Components

for i=1:2:length(varargin)
    switch varargin{i}
        case 'Uis'
            Uis = varargin{i+1}; % m/s
        case 'Uis_bias'
            Uis_bias = varargin{i+1}; %m/s
        case 'P1'
            P1 = varargin{i+1}; % torr
        case 'P1_bias'
            P1_bias = varargin{i+1}; %torr
        case 'T1'
            T1 = varargin{i+1}; % K
        case 'T1_bias'
            T1_bias = varargin{i+1}; %K
        case 'testGasSpec'
            testGasSpec = varargin{i+1}; % cell array of species in test gas
        case 'XtestGas'
            XtestGas = varargin{i+1}; % vector of test gas mole fractions
        case 'XtestGas_bias'
            XtestGas_bias = varargin{i+1}; %bias for fuel mole fraction
        case 'solnMethod'
            solnMethod = varargin{i+1}; % specifies EE, FE, or FF solution procedure
        case 'nasaFile'
            nasaFile = varargin{i+1}; % string name of *.dat file containing NASA polynomial coefficients
        otherwise
            % do nothing
    end
end

[Tout, Pout, ~, ~, ~, ~, ~, ~] = frosh('Uis',Uis,'T1',T1,'P1',P1,'testGasSpec',testGasSpec,'XtestGas',XtestGas,'solnMethod',solnMethod,'nasaFile',nasaFile);
T5_baseline = Tout(2); %[K]
T2_baseline = Tout(1); %[K]
P5_baseline = Pout(2)/101325; %[atm]
P2_baseline = Pout(1)/101325; %[atm]

for ii=1:2
    
    if mod(ii,2)==1
        Uis_temp = Uis+Uis_bias;
    else
        Uis_temp = Uis-Uis_bias;
    end
    
    [Tout, Pout, ~, ~, ~, ~, ~, ~] = frosh('Uis',Uis_temp,'T1',T1,'P1',P1,'testGasSpec',testGasSpec,'XtestGas',XtestGas,'solnMethod',solnMethod,'nasaFile',nasaFile);
    T5_Uis(ii) = Tout(2); %[K]
    T2_Uis(ii) = Tout(1); %[K]
    P5_Uis(ii) = Pout(2)/101325; %[atm]
    P2_Uis(ii) = Pout(1)/101325; %[atm]
    
end

T5_Uis_bias = max(abs(T5_Uis-T5_baseline));
T2_Uis_bias = max(abs(T2_Uis-T2_baseline));
P5_Uis_bias = max(abs(P5_Uis-P5_baseline));
P2_Uis_bias = max(abs(P2_Uis-P2_baseline));


for ii=1:2
    
    if mod(ii,2)==1
        T1_temp = T1+T1_bias;
    else
        T1_temp = T1-T1_bias;
    end
    
    [Tout, Pout, ~, ~, ~, ~, ~, ~] = frosh('Uis',Uis,'T1',T1_temp,'P1',P1,'testGasSpec',testGasSpec,'XtestGas',XtestGas,'solnMethod',solnMethod,'nasaFile',nasaFile);
    T5_T1(ii) = Tout(2); %[K]
    T2_T1(ii) = Tout(1); %[K]
    P5_T1(ii) = Pout(2)/101325; %[atm]
    P2_T1(ii) = Pout(1)/101325; %[atm]
end

T5_T1_bias = max(abs(T5_T1-T5_baseline));
T2_T1_bias = max(abs(T2_T1-T2_baseline));
P5_T1_bias = max(abs(P5_T1-P5_baseline));
P2_T1_bias = max(abs(P2_T1-P2_baseline));

for ii=1:2
    
    if mod(ii,2)==1
        P1_temp = P1+P1_bias;
    else
        P1_temp = P1-P1_bias;
    end
    
    [Tout, Pout, ~, ~, ~, ~, ~, ~] = frosh('Uis',Uis,'T1',T1,'P1',P1_temp,'testGasSpec',testGasSpec,'XtestGas',XtestGas,'solnMethod',solnMethod,'nasaFile',nasaFile);
    T5_P1(ii) = Tout(2); %[K]
    T2_P1(ii) = Tout(1); %[K]
    P5_P1(ii) = Pout(2)/101325; %[atm]
    P2_P1(ii) = Pout(1)/101325; %[atm]
    
end

T5_P1_bias = max(abs(T5_P1-T5_baseline));
P5_P1_bias = max(abs(P5_P1-P5_baseline));
T2_P1_bias = max(abs(T2_P1-T2_baseline));
P2_P1_bias = max(abs(P2_P1-P5_baseline));
%for ii=1:2
%    
%    if mod(ii,2)==1
%        XtestGas_temp = [XtestGas(1)+XtestGas_bias XtestGas(2) XtestGas(3)-XtestGas_bias];
%        %XtestGas_temp = [XtestGas(1)+XtestGas_bias XtestGas(2)-XtestGas_bias];
%    else
%        XtestGas_temp = [XtestGas(1)-XtestGas_bias XtestGas(2) XtestGas(3)+XtestGas_bias];
%        %XtestGas_temp = [XtestGas(1)-XtestGas_bias XtestGas(2)+XtestGas_bias];
%    end
%    
%    [Tout, Pout, ~, ~, ~, ~, ~, ~] = frosh('Uis',Uis,'T1',T1,'P1',P1,'testGasSpec',testGasSpec,'XtestGas',XtestGas_temp,'solnMethod',solnMethod,'nasaFile',nasaFile);
%    T5_XtestGas(ii) = Tout(2);
%    P5_XtestGas(ii) = Pout(2)/101325;
%    
%end


%This is still wrong.  You need to do the full bias analysis.
Num_Components=length(testGasSpec);

%Assumes last gas added had greatest bias, and all other gas biases are a
%calculated as a percent of the greatest bias. Still don't think this is
%right, but probably more conservative than previous method.
for ii=1:2
    if mod(ii,2)==1
switch Num_Components
    case 1 
        XtestGas_temp=[XtestGas(1)];
    case 2 
        XtestGas_temp = [XtestGas(1)+XtestGas_bias(2) XtestGas(2)-XtestGas_bias(2)];
    case 3 
        XtestGas_temp = [XtestGas(1)+((XtestGas_bias(3)*XtestGas(1))/sum(XtestGas(1:2))) XtestGas(2)+((XtestGas_bias(3)*XtestGas(2))/sum(XtestGas(1:2))) XtestGas(3)-XtestGas_bias(3)];
    case 4
        XtestGas_temp = [XtestGas(1)+((XtestGas_bias(4)*XtestGas(1))/sum(XtestGas(1:3))) XtestGas(2)+((XtestGas_bias(4)*XtestGas(2))/sum(XtestGas(1:3))) XtestGas(3)+((XtestGas_bias(4)*XtestGas(3))/sum(XtestGas(1:3))) XtestGas(4)-XtestGas_bias(4)];
    otherwise
            % do nothing
end
    else
switch Num_Components
    case 1
        XtestGas_temp=[XtestGas(1)];
    case 2
        XtestGas_temp = [XtestGas(1)-XtestGas_bias(2) XtestGas(2)+XtestGas_bias(2)];
    case 3
        XtestGas_temp = [XtestGas(1)-((XtestGas_bias(3)*XtestGas(1))/sum(XtestGas(1:2))) XtestGas(2)-((XtestGas_bias(3)*XtestGas(2))/sum(XtestGas(1:2))) XtestGas(3)+XtestGas_bias(3)];
    case 4
        XtestGas_temp = [XtestGas(1)-((XtestGas_bias(4)*XtestGas(1))/sum(XtestGas(1:3))) XtestGas(2)-((XtestGas_bias(4)*XtestGas(2))/sum(XtestGas(1:3))) XtestGas(3)-((XtestGas_bias(4)*XtestGas(3))/sum(XtestGas(1:3))) XtestGas(4)+XtestGas_bias(4)];
    otherwise
            % do nothing
end
    end
    [Tout, Pout, ~, ~, ~, ~, ~, ~] = frosh('Uis',Uis,'T1',T1,'P1',P1,'testGasSpec',testGasSpec,'XtestGas',XtestGas_temp,'solnMethod',solnMethod,'nasaFile',nasaFile);
    T5_XtestGas(ii) = Tout(2); %[K]
    T2_XtestGas(ii) = Tout(1); %[K]
    P5_XtestGas(ii) = Pout(2)/101325; %[atm]
    P2_XtestGas(ii) = Pout(1)/101325; %[atm]
   
end

T5_XtestGas_bias = max(abs(T5_XtestGas-T5_baseline));
T2_XtestGas_bias = max(abs(T2_XtestGas-T2_baseline));
P5_XtestGas_bias = max(abs(P5_XtestGas-P5_baseline));

[Tout, Pout, ~, ~, ~, ~, ~, ~] = frosh_for_sens('Uis',Uis,'T1',T1,'P1',P1,'testGasSpec',testGasSpec,'XtestGas',XtestGas,'solnMethod',solnMethod,'nasaFile',nasaFile);
T5_NASA_bias = Tout(2)-T5_baseline; %[K]
T2_NASA_bias = Tout(1)-T2_baseline; %[K]
P5_NASA_bias = Pout(2)/101325-P5_baseline; %[atm]
P2_NASA_bias = Pout(1)/101325-P2_baseline; %[atm]

B_T5 = sqrt((T5_Uis_bias/T5_baseline)^2 + (T5_T1_bias/T5_baseline)^2 + (T5_P1_bias/T5_baseline)^2 + (T5_XtestGas_bias/T5_baseline)^2 + (T5_NASA_bias/T5_baseline)^2);
B_P5 = sqrt((P5_Uis_bias/P5_baseline)^2 + (P5_T1_bias/P5_baseline)^2 + (P5_P1_bias/P5_baseline)^2 + (P5_XtestGas_bias/P5_baseline)^2 + (P5_NASA_bias/P5_baseline)^2);
B_T2 = sqrt((T2_Uis_bias/T2_baseline)^2 + (T2_T1_bias/T2_baseline)^2 + (T2_P1_bias/T2_baseline)^2 + (T2_XtestGas_bias/T2_baseline)^2 + (T2_NASA_bias/T2_baseline)^2);
B_P2 = sqrt((P2_Uis_bias/P2_baseline)^2 + (P2_T1_bias/P2_baseline)^2 + (P2_P1_bias/P2_baseline)^2 + (P2_XtestGas_bias/P2_baseline)^2 + (P2_NASA_bias/P2_baseline)^2);


T5 = T5_baseline;
P5 = P5_baseline;

end