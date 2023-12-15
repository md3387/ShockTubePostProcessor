%Parse Element name strings to find the names and quantities of individual elements, 
%Handles chemical names w/o coefficients behind each element (Like in CH4 or
%N2O
%M. Hageman 11/2018
%Primary code taken verbatim from Fangjun Jiang on MATLAB Answers:
%https://www.mathworks.com/matlabcentral/answers/13600-how-to-extract-info-from-a-chemical-formula
%Accessed 11/6/2018

%Inputs: 
    %str = Chemical name
%Output: 
    %EleList = cell array of element names
    %NumList = array of element coefficients
    function [EleList NumList]=ParseElementString(str)
%str='C22H10PuCrN2O5'

%remainder of code is verbatim from Fanjun Jiang.
[EleList,Trash,EleEnd]=regexp(str,['[','A':'Z','][','a':'z',']?'],'match');
[Num,NumStart]=regexp(str,'\d+','match');
NumList=ones(size(EleList));
Index=ismember(EleEnd+1,NumStart);
NumList(Index)=cellfun(@str2num,Num);