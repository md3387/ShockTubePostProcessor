function[PlotData, PointData, PlotLables, PointLables]=MlappLasTrace(ShockFilePath,VacuumFilePath,NumHeaderLines,timeColumn,PitchColumn,CatchColumn)
% "MlappLASTrace" - Mitchell D. Hageman October 2024
% PURPOSE:
%   *1) Convert LAS detector voltages into absorbance (absorbance=-ln(It/I0)) and plot.
%   *2) Report average absorbance values at States 1, 2, and 5
%   *3) Use schlieren spikes to identify time [ms] of passage of incident and reflected shocks
%   *4) Shock passage determines test section optical plane State:
%   * State 1 = pre-shock. initial signals seen at beginning of trace.
%   * State 2 = post-incident-shock. Between incident shock wave and contact surface.  SS signal after first schlieren spike.
%   * State 3 = post-rarefaction wave gas.  Not relevant here.
%   * State 4 = pre-shock driver gas.  Not relevant here.
%   * State 5 = post-reflected-shock. Between reflected shock wave and endwall.  SS pressure after Second pressure spike.
% PREREQUISITES
%   * MlappLASTraceHelp.pdf is needed in the PWD (present working directory) if you want the help button to work.
%   * MATLAB Signal Processing Toolbox is needed for filter design.
% INPUTS:
%   * VacuumFilePath - full file path and name of csv or excel file with your Vacuum sample data
%       -example:  'C:\Users\mitchell.hageman\Desktop\Data\20240923_001_VacuumData.csv'
%       -assumes that any dark (zero) signal offset correction has already been applied to the recorded voltages.
%   * ShockFilePath - full file path and name of csv or excel file with your sample data
%       -example:  'C:\Users\mitchell.hageman\Desktop\Data\20240923_001_ShockData.csv'
%       -assumes that any dark (zero) signal offset correction has already been applied to the recorded voltages.
%   * NumHeaderLines - number of header lines in csv file before data begins.
%       -example: My data has two header lines.  Row 1 is data labels, and Row 2 is the offset voltage. So my data  begins in row 3.
%   *timeColumn - column of csv where time trace is.  Mine is Column A so Timecolumn=1.
%   *PitchColumn - column of csv where reference detector trace is.  Mine is column C so PitchColumn=3.
%   *CatchColumn - column of csv where signal detector  trace is.  Mine is column D so CatchColumn=4.

% OUTPUTS:
%  * MatrixName 
%       -Variable name [units]                      -DataType   Notes
%----------------------------------------------------------------------------------------------------------
%   * PlotData
%       - time [ms],                                -Vector
%       - absorbance [-],                           -Vector    = -ln(It/IoAvg)
%       - Io [-],                                   -Vector    = VacuumSignalVoltage./VacuumReferenceVoltage
%       - It [-],                                   -Vector    = SampleSignalVoltage./SampleReferenceVoltage
%       - SampleReferenceVoltage [V]                -Vector  
%       - Sample SignalVoltage [V]                  -Vector
%   *PointData
%       -State 1 absorbance average[-]             - Value
%       -State 1 absorbance standard deviation [-] - Value
%       -State 2 absorbance average[-]             - Value
%       -State 2 absorbance standard deviation [-] - Value
%       -State 5 absorbance average[-]             - Value
%       -State 5 absorbance standard deviation [-] - Value
%       -t_incident [ms]                           - Value  = time incident shock passed optical plane
%       -t_reflected [ms]                          - Value  = time reflected shock passed optical plane
%   * PlotLables - Column headers for PlotData     - Character Matrix
%   * PointLables - Data lables for PointData      - Character Matrix
%   *"Write Plot Data to File" Button writes PlotData and Point Data to .csv in user-selected folder, with lables in the first row.
% VERSION NUMBER:
%   * 1.0: October 2024 - initial release, Mitchell D. Hageman
%% Read in detector voltages
ShockData = readmatrix(ShockFilePath,'NumHeaderLines',NumHeaderLines);
time = ShockData(:,timeColumn); %[ms]
SampleReferenceVoltage=ShockData(:,PitchColumn); %[V]
SampleSignalVoltage=ShockData(:,CatchColumn); %[V]

VacuumData= readmatrix(VacuumFilePath,'NumHeaderLines',NumHeaderLines);
VacuumReferenceVoltage=VacuumData(:,PitchColumn); %[v]
VacuumSignalVoltage=VacuumData(:,CatchColumn); %[v]

It=SampleSignalVoltage./SampleReferenceVoltage; %It -with CMR.(Assumes Dark signal offset is already applied to loaded data)
Io=VacuumSignalVoltage./VacuumReferenceVoltage;%Io -with CMR.(Assumes Dark signal offset is already applied to loaded data)
IoAvg=mean(Io);
IoStdev=std(Io);
absorbance = real(-log((It)./IoAvg));% With Common mode rejection (Assumes Dark signal offset is already applied to loaded data)


%% Create Interactive Figure and ID State Boundaries
fig = uifigure('Name', 'Interactive Plot', 'Position', [500 300 800 600]); % Create the main GUI window (figure)

% Create a button that opens a help file
helpButton = uibutton(fig, 'push', 'Text', 'Help!', ...
    'Position', [600 10 150 25], ...
    'ButtonPushedFcn', @(btn, event) openHelpFile());

% Create axes for plotting
ax = uiaxes(fig, 'Position', [50 100 700 475]);
plot(ax,time,SampleReferenceVoltage,'DisplayName',"Reference Voltage [V]");
hold(ax, 'on');
plot(ax,time,SampleSignalVoltage,'DisplayName',"Signal Voltage [V]");
legend(ax, 'show');
ax.YLabel.String ="Voltage [V]";
ax.XLabel.String = 'time [ms]';
hold(ax, 'off');
proceedButton = uibutton(fig, 'push', 'Text', 'Proceed to It and Io plot?', 'Position', [200 10 250 25], ...
    'ButtonPushedFcn', @(btn, event) resumeCalculation());
% Pause execution here until the user presses the Proceed button
uiwait(fig); % Pauses the code execution

cla(ax);
plot(ax,time,Io,'DisplayName',"Io [-]");
hold(ax, 'on');
plot(ax,time,It,'DisplayName',"It [-]");
legend(ax, 'show');
ax.YLabel.String ="Intensity [-]";
ax.XLabel.String = 'time [ms]';
hold(ax, 'off');
proceedButton = uibutton(fig, 'push', 'Text', 'Proceed to absorbance plot?', 'Position', [200 10 250 25], ...
    'ButtonPushedFcn', @(btn, event) resumeCalculation());
% Pause execution here until the user presses the Proceed button
uiwait(fig); % Pauses the code execution

cla(ax);
plot(ax,time,absorbance,'DisplayName',"absorbance [-]"); %ax=axes you created, time=x, Pressure_atm=y, DisplayName=legend command, "Pressure [atm]" = display name.
hold(ax, 'on');
legend(ax, 'show');
ax.YLabel.String ="Absorbance [-]";
ax.XLabel.String = 'time [ms]';

State1EndLabel = uilabel(fig, 'Text', 'State 1 End [ms]:', 'Position', [20 75 100 25]); %Create Label
State1EndInput = uieditfield(fig, 'numeric', 'Position', [120 75 25 25], 'Value', 0); %Create input box

State2StartLabel = uilabel(fig, 'Text', 'State 2 Start [ms]:', 'Position', [150 75 100 25]); %Create Label
State2StartInput = uieditfield(fig, 'numeric', 'Position', [250 75 25 25], 'Value', 0); %Create input box
State2EndLabel = uilabel(fig, 'Text', 'State 2 End [ms]:', 'Position', [280 75 100 25]); %Create Label
State2EndInput = uieditfield(fig, 'numeric', 'Position', [380 75 25 25], 'Value', 0); %Create input box

State5StartLabel = uilabel(fig, 'Text', '.State 5 Start [ms]', 'Position', [410 75 100 25]); %Create Label
State5StartInput = uieditfield(fig, 'numeric', 'Position', [510 75 25 25], 'Value', 0); %Create input box
State5EndLabel = uilabel(fig, 'Text', 'State 5 End [ms]', 'Position', [540 75 100 25]); %Create Label
State5EndInput = uieditfield(fig, 'numeric', 'Position', [640 75 25 25], 'Value', 0);

% Create a button that will allow the user to proceed after entering the value
proceedButton = uibutton(fig, 'push', 'Text', 'Calculate Averages and Shock Passage', 'Position', [200 10 250 25], ...
    'ButtonPushedFcn', @(btn, event) resumeCalculation());
% Pause execution here until the user presses the Proceed button
uiwait(fig); % Pauses the code execution

%Convert User Input Times to Indices
State1End=State1EndInput.Value; %Record Value [ms]
[~, State1EndIndex]=min(abs(time-State1End)); %[index]
State2Start=State2StartInput.Value; %Record Value [ms]
[~, State2StartIndex]=min(abs(time-State2Start)); %[index]
State2End=State2EndInput.Value; %Record Value
[~, State2EndIndex]=min(abs(time-State2End)); %[index]
State5Start=State5StartInput.Value; %Record Value
[~, State5StartIndex]=min(abs(time-State5Start)); %[index]
State5End=State5EndInput.Value; %Record Value
[~, State5EndIndex]=min(abs(time-State5End)); %[index]

%Find Average Measured and Derived Quantities between the Indices
State1It=mean(It(1:State1EndIndex));
State1ItStdev=std(It(1:State1EndIndex));
State1Absorbance=mean(absorbance(1:State1EndIndex));
State1AbsorbanceStdev=std(absorbance(1:State1EndIndex));

State2It=mean(It(State2StartIndex:State2EndIndex));
State2ItStdev=std(It(State2StartIndex:State2EndIndex));
State2Absorbance=mean(absorbance(State2StartIndex:State2EndIndex));
State2AbsorbanceStdev=std(absorbance(State2StartIndex:State2EndIndex));

State5It=mean(It(State5StartIndex:State5EndIndex));
State5ItStdev=std(It(State5StartIndex:State5EndIndex));
State5Absorbance=mean(absorbance(State5StartIndex:State5EndIndex));
State5AbsorbanceStdev=std(absorbance(State5StartIndex:State5EndIndex));

%Define t_incident and t_reflected based on maximum gradients between state definitions
[~,t_IncidentIndex] = max(abs(absorbance(State1EndIndex:State2StartIndex)));
t_incident = time(t_IncidentIndex+State1EndIndex);%[ms] time of reflected shock arrival at sidewall.
plot(ax,[t_incident t_incident],[0 max(absorbance)],'DisplayName','t_{incident}');

[~,t_ReflectedIndex] = max(abs(absorbance(State2EndIndex:State5StartIndex)));
t_reflected = time(t_ReflectedIndex+State2EndIndex);%[ms] time of reflected shock arrival at sidewall.
plot(ax,[t_reflected t_reflected],[0 max(absorbance)],'DisplayName','t_{reflected}');


PlotData=[time, absorbance, Io, It, SampleReferenceVoltage, SampleSignalVoltage];
PlotLables={'time [ms]', 'absorbance [-]', 'Io [-]', 'It [-]','Pitch [V]', 'Catch [V]'};
PointData=[IoAvg, IoStdev, State1It, State1ItStdev, State1Absorbance, State1AbsorbanceStdev,...
                           State2It, State2ItStdev, State2Absorbance, State2AbsorbanceStdev,...
                           State5It, State5ItStdev, State5Absorbance, State5AbsorbanceStdev...
                           t_incident, t_reflected];
PointLables={'IoAvg', 'Io Stdev','State1 It', 'State1 It Stdev', 'State1 Absorbance', 'State1 AbsorbanceStdev',...
                                'State2 It', 'State2 It Stdev', 'State2 Absorbance', 'State2 Absorbance Stdev',...
                                'State5 It', 'State5 It Stdev', 'State5 Absorbance', 'State5 Absorbance Stdev',...
                                't_incident[ms]', 't_refelcted [ms]'};
%Create a button that will save data to csv
plotButton = uibutton(fig, 'push', 'Text', 'Write Data to CSV', ...
    'Position', [20 10 150 25], ...
    'ButtonPushedFcn', @(btn, event) WritePlotDataToCSV(PlotData,PlotLables,PointData,PointLables),...
    'BackgroundColor', [0, 1, 0]);% RGB color for button background;


%% Callback Functions

%Callback Function to write data
    function WritePlotDataToCSV(PlotData,PlotLables,PointData,PointLables)
        [PlotDatafile,location]= uiputfile('.xlsx');
        writematrix(PlotData,[char(location),char(PlotDatafile)],'Range','A2');
        writecell(PlotLables,[char(location),char(PlotDatafile)],'Range','A1');
        writematrix(PointData,[char(location),char(PlotDatafile)],'Range','G2');
        writecell(PointLables,[char(location),char(PlotDatafile)],'Range','G1');
    end

%Callback function to open the help file
    function openHelpFile()
        helpFile = 'MlappLASTraceHelp.pdf';  % Path to the help file
        if isfile(helpFile) % Use MATLAB's open function to open the file
            open(helpFile);
        else
            % Display a message if the file doesn't exist
            uialert(fig, 'Help file not found!', 'Error');
        end
    end

%Callback function for "Proceed..." button
    function resumeCalculation()
        % Resume program execution
        uiresume(fig);
    end

hold(ax, 'off'); %allow next plug-in to move on to the next figure.
end