function[BifurcationHeight]=PredictBifurcation(Mout, gammaOut, MW_mix)
%Ref 1: THE INTERACTION OF A REFLECTED SHOCK WAVE WITH THE BOUNDARY LAYER IN A SHOCK TUBE 
%       Herman Mark (NACA TM 1418) 1958
%Ref 2: Measurement of reflected-shock bifurcation over a wide range of gas composition and pressure
%       E. L. Petersen · R. K. Hanson, Shock Waves (2006) 15:333–340DOI 10.1007/s00193-006-0032-3

M_3=sqrt((2*gammaOut(1)*Mout(1)^2-(gammaOut(1)-1))/((gammaOut(1)-1)*Mout(1)^2+2));%[-] Mach number at State 3 (Ref 1, Eq II-23)
M_BL=(2*(gammaOut(1)-1)*Mout(1)^2+(3-gammaOut(1)))/((gammaOut(1)+1)*Mout(1));%[-] Mach number in the LAMINAR BOUNDARY LAYER (Ref 1, Eq III-3)
PBLoverP2=(1+((gammaOut(1)-1)/2)*M_BL^2)^(gammaOut(1)/(gammaOut(1)-1));%[-] Ref1 eq III-15
P5overP2=((2*gammaOut(1)/(gammaOut(1)+1))*M_3^2)-((gammaOut(1)-1)/(gammaOut(1)+1)); %[-] Ref1 eq III-14
Mrs=Mout(2);%Reflected shock Mach number
if PBLoverP2<=P5overP2
   BifurcationHeight=7.5*(Mout(1)^1.07)/(gammaOut(1)^2.66*MW_mix^0.37); % [mm] Ref2 eq 10 
   % "The thickness of the boundary layer...was neglected... this extra height should be added.
   
   %{
   %See Ref 2 Figure 1 for definition of terms 
   P_OBL=? %Pressure in the boundary layer at O?
   t_O %time where bifurcation starts (Determine from 1st pressure rise)
   t_A %time where shock passes. (Determine from Schlieren spike?)
   deltatao=t_sp-t_bs;
   x1 = deltatao*Urs; %[mm] Ref2 eq 4
   SinSquared_theta_1=((gammaOut(?)+1)*(P_OBL/Pout(1))+gamma(?)-1)/(2*gamma(?)*Mrs^2)% Ref2 eq7 rearranged.
   BifurcationHeight2=x1*tan(theta_1); %[degrees] Ref2 eq5
   BifurcationHeightUncertainty=(BifurcationHeight-BifurcationHeight2);% [mm]
   delta= %Need "Influence of Reflected Shock and Boundary‐Layer Interaction on Shock‐Tube Flows" L. Davies; J. L. Wilson Phys. Fluids 12, I-37–I-43 (1969)
   alpha= %Need "Influence of Reflected Shock and Boundary‐Layer Interaction on Shock‐Tube Flows" L. Davies; J. L. Wilson Phys. Fluids 12, I-37–I-43 (1969)
   theta_2=???
   %}

else
    BifurcationHeight='NA';
end
