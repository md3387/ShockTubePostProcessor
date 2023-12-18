# ShockTubePostProcessor
Post-Processor for a combustion shock tube

File naming and structure:
---------------------------------------------------------------------------------------------------------------------------------------------------------------------
-All files are .csv files
-All files are named in the same way: date(yyyymmdd) _ run#(00#) _ file type.
-Date ex: 20230929 is September 29, 2023)
-Run# ex: run one is _001, run ten is _010)

The same 4 file types are output for every experiment:
yyyymmdd_run_Experiment Parameters
yyyymmdd_run_DarkFullData
yyyymmdd_run_TimerData
yyyymmdd_run_ShockData

***The content of these file types is:***

yyyymmdd_run_Experiment Parameters: 
1KB, Describes tube initial conditions. 5 rows, two columns
Ex contents: 
Vacuum Pressure (torr),	3.80E-06
Test Gas Pressure (torr),	3.315
Region 1 Pressure (torr),	3.315
Region 1 Temperature (degC),	25.6
Region 4 Pressure (psi),	7.2

yyyymmdd_run_DarkFullData
~8,100KB, Full Oscilloscope record for static tube under vacuum. It's used for offset-correction, so all nearly all data streams look like noise around a near-zero signal
-100,002 rows, 9 columns. Row 1 is Headers, Row 2 is an offset, next 100,000 rows are data
- Column 1 is time (ms), Columns 2-9 are Oscope channels 0-7
Headers:
Time (ms)	Kistler Pressure	OH Side	HeNe Pitch	HeNe Catch	OH End	PCB 5	CO2 Pitch	CO2 Catch


yyyymmdd_run_TimerData
1KB, Counter-Timer data used to determine shock velocity.
5 rows, 4 columns (but rows are likely to be added if additional counter-timers are added)
Ex contents:
Channel Name	 Start Source	 End Source	 Measured Time (s)
Counter 0	 PCB 1	 PCB 2	 2.683700E-4
Counter 1	 PCB 2	 PCB 3	 2.668700E-4
Counter 2	 PCB 3	 PCB 4	 2.674800E-4
Counter 3	 PCB 4	 PCB 5	 2.708100E-4

yyyymmdd_run_ShockData
~8,100KB, Full Oscilloscope record for experiment. This is the primary experimental data, and all the interesting O-scope traces are found here. 
Identical structure to yyyymmdd_run_DarkFullData




I've included a few days of experiment files as example data. This is not an exhaustive list and does not necessarily represent every type of experiment that could be experienced.  However, it is a good start to understanding the format of the Oscope traces.

Example Runs and details:
-----------------------------------------------------------------------------------------------------------------------------------------------------------------------
yyyymmdd_run# -  Trigger - Ignition  -   Notes                                                        Features to identify from "Kistler Pressure" trace:
20230920_001  -    1ms   -    Fast   -   dpdt does not have peaks near incident or reflected shock,   1)incident shock, 2)Reflected shock 3)shock bifurcation 4)ignition
20230920_001  -    1ms   -  Medium   -   dpdt has weak peaks near incident and reflected shocks,      1)incident shock, 2)Reflected shock 3)shock bifurcation 4)ignition
20230929_001  -    2ms   - very Fast -   really weak incident shock, non-obvious shock bifurcation,   1)incident shock, 2)Reflected shock 3)shock bifurcation 4)ignition
20230929_005  -    2ms   -   None    -   weak incident shock, no shock bifurcation,                   1)incident shock, 2)Reflected shock  
20231116_006  -    2ms   -   None    -   but reflected shock is offset, no shock bifurcation,         1)incident shock, 2)Reflected shock  
