%Use parsed CHON molecule elements to find MW
%M. Hageman 11/2018

%Inputs:
%SpeciesName example: 'CH4' or 'O2' etc.  Only CHON+Ar+He containing molecules will work.

%Outputs:
%MW - molecular weight of the input species [kg/kmol or g/mol]

%Requires Subfunctions:
%ParseElementString(str)

function [MW] = CHON_MW(SpeciesName)

[ElementList, NumberList]=ParseElementString(SpeciesName);
CHON = {'C' 'H' 'O' 'N' 'Ar' 'He'}; %strings for the four elements we're concerned about plus two inerts.
CHON_MW = [12.011 1.008 15.999 14.007 39.948 4.003]; %molecular weights for respective eleements [kg/kmol] or [g/mol]

for i=1:length(NumberList) %i = number of elements in the molecule string
    [~,loc(i)]=ismember(ElementList(i),CHON); %match Element 1 string to a string in the CHON String array. The index of that string matches the index of the MW in the CHON_MW array.
    MW_Element(i)=CHON_MW(loc(i)); %output Element MW from the CHON_MW array
    MW_Atoms(i)=MW_Element(i)*NumberList(i); %multiply element MW by the number of atoms of that element in the molecule
end
MW = sum(MW_Atoms); %[kg/kmol] or [g/mol]



