function[CT_file, CT_data,shock_file, shock_data, dark_file, dark_data, EP_file, EP_data]=GetFiles(file_base, date, run_no)

filestring=strcat(file_base,num2str(date),'_',num2str(run_no,'%03d'),'_');
CT_file = strcat(filestring,'TimerData.csv');
CT_data_full = readmatrix(CT_file,'NumHeaderLines',1);
CT_data = CT_data_full(:,end);

shock_file = strcat(filestring,'ShockData.csv'); 
shock_data = readmatrix(shock_file,'NumHeaderLines',2); %Load and extract shock data

dark_file = strcat(filestring,'DarkFullData.csv'); 
dark_data = readmatrix(dark_file,'NumHeaderLines',2); %Load and extract dark data

EP_file = strcat(filestring,'ExperimentParameters.csv');
EP_data = readmatrix(EP_file,'NumHeaderLines',2);%Load and extract experimental parameters
