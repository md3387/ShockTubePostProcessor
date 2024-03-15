function [EndwalltoMidpoints,  U, xx, Upred, UpredCI, Uis, Uis_bias, time_to_EW, UatPCB, comments]=MlappShockVelocity(CT_data,DMdist,numberOfCounters,firstDMpcbPort,EndwalltoTestSection, comments)

%Preallocate matrices
CtrDist=zeros(1,numberOfCounters);
EndwalltoPCB=zeros(1,numberOfCounters+1);
EndwalltoMidpoints=zeros(1,numberOfCounters);
weight=ones([1,numberOfCounters]); 
UatPCB=zeros(1,numberOfCounters);
U=zeros(1,numberOfCounters); 

for i=1:numberOfCounters
CtrDist(i)=DMdist(i+firstDMpcbPort-1);
EndwalltoPCB(1)=61.9*0.0254;% [m] Distance from PCB to endwall (specified. Ref.2 CST-1005, CST-1006, CST-1007)
EndwalltoPCB(i+1)=EndwalltoPCB(i)-CtrDist(i);% [m] Distance from PCB to endwall
EndwalltoMidpoints(i)=EndwalltoPCB(i)-(CtrDist(i)/2);% [m] Distance from port spacing midpoint to endwall
U(i)=CtrDist(i)/CT_data(i);
end

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
Uis_CI = confint(U_fit);
Uis_bias = (Uis_CI(2,2)-Uis_CI(1,2))/2;
%Weighted nonlinear regression model, 1st degree polynomial, Weights all =1, so should give same result as above, but easier to plot Confidence intervals
[model, S]=polyfit(EndwalltoMidpoints,U,1);
start=[model(1);model(2)];
modelstring=@(b,EndwalltoMidpoints) b(1).*EndwalltoMidpoints+b(2);
U_wnlm=fitnlm(EndwalltoMidpoints,U,modelstring,start,'Weight',weight);
xx=linspace(0,2)';
[Upred,UpredCI]=predict(U_wnlm,xx,'Simultaneous',true);

UatIncident = U_fit.p1*EndwalltoTestSection + U_fit.p2; %[m/s] shock velocity at T_incident?
time_to_EW = roots([0.5*U_fit.p1 UatIncident -EndwalltoTestSection]); %[s]
time_to_EW = time_to_EW(2)*1000; %convert from [s] to [ms]


%Calculate time of incident shock arrival at PCBs
for j=1:numberOfCounters+1
UatPCB(i) = U_fit.p1*EndwalltoPCB(i) + U_fit.p2;
end

%{
figureIndex=figureIndex+1;
figure(figureIndex)
hold on
plot(EndwalltoMidpoints(1),U(1),'r O')
plot(EndwalltoMidpoints(2),U(2),'r O')
plot(EndwalltoMidpoints(3),U(3),'r O')
plot(EndwalltoMidpoints(4),U(4),'r O')
plot(U_fit)
%%line(xx,predict(U_wnlm,xx),'color','b','linestyle', '--')
plot(xx,Upred,'b--',xx, UpredCI, 'b--','DisplayName','WNLM Fit & CI')
plot(0, Uis, 'b O', 'MarkerFaceColor', 'blue', 'DisplayName','Uis Endwall')
plot(0, Uis+Uis_bias, 'b O', 'DisplayName','Uis Bias')
plot(0, Uis-Uis_bias, 'b O', 'DisplayName','Uis Bias')
xlabel('Distance From Endwall [m]')
ylabel('Shock Velocity [m/s]')
hold off
%}
