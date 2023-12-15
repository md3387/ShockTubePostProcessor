%Mole fractions of multiple equivalence ratios

%Note: Added MolarFraction(i) on 7/11/2019, so if dependent functions stop
%working later, that may be the issue.

%Inputs
%Species = {'CH4' 'O2' 'N2'};
%PhiRange = [0.7 0.1 1.5]; [low step high]
%MolesInertperMoleOxidizer = 3.76; for air

%Uses Functions
%CHON_MW.m
%ParseElementString.m
%FindStoich.m

function [MoleFractionArray] = PhiToMoles(Species, MolesInertperMoleOxidizer, PhiRange)


[OF_Stoich, MolArrayStoich] = FindStoich(Species, MolesInertperMoleOxidizer)
i=0
for Phi=PhiRange(1):PhiRange(2):PhiRange(3)
    i=i+1;
    ReportPhi(i)=Phi;
    FuelMoles(i)=1;
    OxidizerMoles(i)=MolArrayStoich(2)/Phi;
    InertMoles(i)=OxidizerMoles(i)*MolesInertperMoleOxidizer;
    MolarFraction(i)=FuelMoles(i)/(OxidizerMoles(i)+InertMoles(i));
end
Phi'

MoleFractionArray=[ReportPhi' FuelMoles' OxidizerMoles' InertMoles' MolarFraction'];