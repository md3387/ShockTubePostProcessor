function[frequencystring,frequencydouble]=FrequencySwitchCase(value)
switch value
    case 'Kistler PXD'
        frequencystring ='100000';
    case 'OH SW'
        frequencystring ='2600';
    case 'HeNe Pitch'
        frequencystring ='1000000';%[Hz] PDAVJ5 - DC-1MHz
    case 'HeNe Catch'
        frequencystring ='1000000';%[Hz] PDAVJ5 - DC-1MHz
    case 'OH EW'
        frequencystring ='100000'; %PMTFrequency=80000000; %[Hz] end wall
    case 'PCB 5'
        frequencystring ='100000';
    case 'CO2 Pitch'
        frequencystring ='1000000';
    case 'CO2 Catch'
        frequencystring ='1000000';
    otherwise
        frequencystring ='';
end
frequencydouble=str2double(frequencystring); %