function [Tout, Pout, rhoOut, gammaOut, VsoundOut, Mout, Urs, phi] = frosh(varargin)
% "FROzen SHock solver"
% PURPOSE: 
%   * Calculate conditions behind incident and reflected shocks in the
%       a conventional shock tube or in a high pressure shock tube
% INPUTS:
%   * Uis = incident shock velocity at endwall (m/s)
%   * T1 = temerature in region 1 (K)
%   * P1 = pressure in region 1 (torr)
%   * testGasSpec = cell array of strings of gases in test gas mixture
%         e.g.: {'AR','O2','C3H8'}
%   * XtestGas = vector of mole fractions of gases in test gas mixture
%         e.g.: [0.96,0.03,0.01]
%   * solnMethod = string specifying vibrational equilibration preference
%        * 'EE' = vibrationally equilibrated in R2 and R5
%        * 'FE' = vibrationally frozen in R2, vibrationally equilibrated in R5
%        * 'FF' = vibrationally frozen in R2 and R5
%   * nasaFile = string name of *.dat file containing NASA polynomial 
%         coefficients for molecules in shock tube gas mixture.  
%         Should be in the same folder as "aerofrosh.m" function 
%         e.g.: 'sandia.dat'
% OUTPUTS:
%   * Tout = [T2 T5]
%        * T2 = temperature in region 2 (K)
%        * T5 = temperature in region 5 (K)  
%   * Pout = [P2 P5]
%        * P2 = pressure in region 2 (Pa)   
%        * P5 = pressure in region 5 (Pa)
%   * rhoOut = [rho1 rho2 rho5]
%        * rho1 = density in region 1 (kg/m3)
%        * rho2 = density in region 2 (kg/m3)
%        * rho5 = density in region 5 (kg/m3)
%   * gammaOut = [gamma1 gamma2 gamma5]
%        * gamma1 = cp1/cv1 in region 1
%        * gamma2 = cp2/cv2 in region 2
%        * gamma5 = cp5/cv5 in region 5
%   * VsoundOut = [Vsound1 Vsound2 Vsound5]
%        * Vsound1 = speed of sound in region 1 (m/s)
%        * Vsound2 = speed of sound in region 2 (m/s)
%        * Vsound5 = speed of sound in region 5 (m/s)
%   * Mout = [Mis Mrs]
%        * Mis = Mach number of incident shock
%        * Mrs = Mach number of reflected shock
%   * Urs = velocity of reflected shock at endwall (m/s)
%   * phi = equivalence ratio of mixture
% DEVELOPMENT HISTORY: 
%   * Originally written by Dr. David Davidson 
%   * Transcribed into MATLAB code and updated by Matthew Campbell.
% VERSION NUMBER:
%   * 1.0: May 2011 - initial release
%   * 1.1: 5/28/2012 - check for convergence
%   * 1.2: 6/2/2012 - better SANDIA.DAT file reading routine

%% PARSE INPUT 
    
% where are we?
    %fprintf('Parsing Input\n')
% Grab information from "varargin"
    for i=1:2:length(varargin)
        switch varargin{i}
            case 'Uis'
                Uis = varargin{i+1}; % m/s
            case 'P1'
                P1 = varargin{i+1}; % torr
            case 'T1'
                T1 = varargin{i+1}; % K
            case 'testGasSpec'
                testGasSpec = varargin{i+1}; % cell array of species in test gas
            case 'XtestGas'
                XtestGas = varargin{i+1}; % vector of test gas mole fractions
            case 'solnMethod'
                solnMethod = varargin{i+1}; % specifies EE, FE, or FF solution procedure
            case 'nasaFile'
                nasaFile = varargin{i+1}; % string name of *.dat file containing NASA polynomial coefficients
            otherwise
                % do nothing
        end
    end
    
%% CHECK INPUTS

% where are we?
    %fprintf('Checking Input\n')
% if vector is column, invert it so that it is a row vector.
    if size(testGasSpec,1)>size(testGasSpec,2)
        testGasSpec=testGasSpec';
    end
    if size(XtestGas,1)>size(XtestGas,2)
        XtestGas=XtestGas';
    end
% check gas vectors
    if (size(testGasSpec,1) ~= size(XtestGas,1)) || (size(testGasSpec,2) ~= size(XtestGas,2))
        error('Test gas species and mole fraction vectors must have same size')
    end
% check temperature
    if T1<273.15
        error('Temperature should be in Kelvin, not Celcius!');
    end
% ensure that mole fractions sum to 1.0
    if abs(sum(XtestGas)-1) > 0.0000000001
        error('Test gas mole fractions do not sum to 1.0');
    end

%% INITIAL COMPUTATIONS

% where are we?
    %fprintf('Performing Initial Computations\n')
% rename and CAPITALIZE species vector
    speciesVec = upper(testGasSpec);
% CAPITALIZE solution method
    solnMethod = upper(solnMethod);
% Transform pressure units
    P1torr = P1; % torr
    P1 = P1torr*(1/760)*(101325/1); % torr --> Pa
% The velocity of the test gas relative to the shock in region 1 (u1) is
% equal to the incident shock velocity (Uis)
    u1 = Uis; % m/s
    
%% CONSTANTS

% where are we?
    %fprintf('Defining Constants\n')
% elemental molecular weights
    eleMolWts = { % g/mol (cell array)
        'AC'	227
        'AG'	107.8682
        'AL'	26.9815386
        'AM'	243
        'AR'	39.948
        'AS'	74.92160
        'AT'	210
        'AU'	196.966569
        'B'     10.811
        'BA'	137.327
        'BE'	9.012182
        'BH'	262
        'BI'	208.98040
        'BK'	247
        'BR'	79.904
        'C'     12.0107
        'CA'	40.078
        'CD'	112.411
        'CE'	140.116
        'CF'	251
        'CL'	35.453
        'CM'	247
        'CO'	58.933195
        'CR'	51.9961
        'CS'	132.9054519
        'CU'	63.546
        'DB'	262
        'DY'	162.500
        'ER'	167.259
        'ES'	252
        'EU'	151.964
        'F'     18.9984032
        'FE'	55.845
        'FM'	257
        'FR'	223
        'GA'	69.723
        'GD'	157.25
        'GE'	72.64
        'H'     1.00794
        'HE'	4.002602
        'HF'	178.49
        'HG'	200.59
        'HO'	164.93032
        'HS'	265
        'I'     126.90447
        'IN'	114.818
        'IR'	192.217
        'K'     39.0983
        'KR'	83.798
        'LA'	138.90547
        'LI'	6.941
        'LR'	262
        'LU'	174.9668
        'MD'	258
        'MG'	24.305
        'MN'	54.938045
        'MO'	95.96
        'MT'	266
        'N'     14.0067
        'NA'	22.98976928
        'NB'	92.90638
        'ND'	144.242
        'NE'	20.1797
        'NI'	58.6934
        'NO'	259
        'NP'	237
        'O'     15.9994
        'OS'	190.23
        'P'     30.973762
        'PA'	231.03588
        'PB'	207.2
        'PD'	106.42
        'PM'	145
        'PO'	209
        'PR'	140.90765
        'PT'	195.084
        'PU'	244
        'RA'	226
        'RB'	85.4678
        'RE'	186.207
        'RF'	261
        'RH'	102.90550
        'RN'	222
        'RU'	101.07
        'S'     32.065
        'SB'	121.760
        'SC'	44.955912
        'SE'	78.96
        'SG'	263
        'SI'	28.0855
        'SM'	150.36
        'SN'	118.710
        'SR'	87.62
        'TA'	180.94788
        'TB'	158.92535
        'TC'	98
        'TE'	127.60
        'TH'	232.03806
        'TI'	47.867
        'TL'	204.3833
        'TM'	168.93421
        'U'     238.02891
        'V'     50.9415
        'W'     183.84
        'XE'	131.293
        'Y'     88.90585
        'YB'	173.054
        'ZN'	65.38
        'ZR'	91.224};
% universal gas constant
    Ru = 8.314471; % J/mol*K
% reference temperature
    Tref = 298.15; % K
    
%% READ SANDIA THERMODYNAMIC DATA FILE, COMPUTE SPECIES MOLECULAR WEIGHTS

% where are we?
    %fprintf('Reading Sandia Thermo Data\n')
% Read file (see subfunction below)
    [coefVec, nameEle, numEle] = readSandia(speciesVec, nasaFile);
% where are we?
    %fprintf('Computing Species Molecular Weights\n')
% Find molecular weight of each species
    specMolWts = getSpecMolecWts(nameEle, numEle, eleMolWts); % g/mol
% where are we?
    %fprintf('Computing Rotational Degrees of Freedom\n')
% Find rotational degrees of freedom for each species
    drotVec = getSpecDrot(speciesVec, nameEle, numEle);

%% DETERMINE GAS PROPERTIES IN REGION ONE

% where are we?
    %fprintf('Determining Gas Properties in Region One\n')
% Weighted sum by the mole fractions in the gas
% multiply then add mole fractions by species molec weights 
    gasMolWt = XtestGas*specMolWts'; % g/mol
% ideal gas constant 
    Rgas = (Ru/gasMolWt)*1000; % (J/mol*K)*(mol/g)*(1000 g/1 kg) = J/kg*K
% specific enthalpy and specific heat of gas in region 1
    if strcmp(solnMethod,'EE')
        % compute from polynomial fits
            [Cp1, h1] = getCph_T(T1, XtestGas, gasMolWt, Ru, coefVec, speciesVec); % J/kg*K and J/kg
    else % 'FE' or 'FF'
        % compute from translational and rotational degrees of freedom
            Cp1moleTR = (Ru + (3/2)*Ru + (1/2)*Ru*drotVec)*XtestGas'; % J/mol*K; Cp = R + Cv; Cv = Trans + Rot; weighted sum by mole fraction
            Cp1TR = (Cp1moleTR/gasMolWt)*1000; % (J/mol*K)*(mol/g)*(1000 g/1 kg) = J/kg*K
        % compute from polynomail fits, including vibration
            [Cp1TRV, h1TRV] = getCph_T(T1, XtestGas, gasMolWt, Ru, coefVec, speciesVec); % J/kg*K and J/kg
        % compute mixture heat of formation (298.15 K)
            [dummy, hForm] = getCph_T(Tref, XtestGas, gasMolWt, Ru, coefVec, speciesVec); % J/kg*K and J/kg
        % compute h1 from specific heat, including only trans and rot
            h1TR = hForm + Cp1TR*(T1-Tref); % J/kg
        % compute vibrational contribution to enthalpy in region 1
        % this is technically an integral of cvdT, which is internal energy rather than enthalpy.  
        % hence the "e" rather than "h" designation
            e1V = h1TRV - h1TR;
        % our effective Cp1 and h1 are what the polynomial fits specified
            Cp1 = Cp1TRV;
            h1 = h1TRV;
    end
% specific volume
    v1 = Rgas*T1/P1; % m3/kg
% density
    rho1 = 1/v1; % kg/m3
% constant volume specific heat
    Cv1 = Cp1-Rgas; %J/kg*K
% gamma
    gamma1 = Cp1/Cv1; % no units
% speed of sound in region one
    Vsound1 = sqrt(gamma1*P1/rho1); % m/s
% mach number of incident shock
    Mis = u1/Vsound1; % no units

%% SOLVE IDEAL INCIDENT SHOCK EQNS

% where are we?
    %fprintf('Solving Incident Shock Conditions\n')
% Before Newton-Rhapson guess using ideal relations:
%     P2/P1 = ((2*gamma1*Mis^2-(gamma1-1))/(gamma1+1))
%     T2/T1 = (gamma1*Mis^2-(gamma1-1)/2)*((gamma1-1)/2*Mis^2+1)/(((gamma1+1)/2)^2*Mis^2)
%	  gamma = cp1/(cp1-1) --> assumes driven gas is mostly Ar with Cp ~ 2.5
%	  Mis = U1/Vsound1
% compute ideal pressure in region 2
    P2ideal = P1*((2*gamma1*Mis^2-(gamma1-1))/(gamma1+1)); % Pa
% compute ideal temperature in region 2
    T2ideal = T1*(gamma1*Mis^2-(gamma1-1)/2)*((gamma1-1)/2*Mis^2+1)/(((gamma1+1)/2)^2*Mis^2); % K
% specific volume
    v2ideal = Rgas*T2ideal/P2ideal; % m3/kg
% Solve the Incident Shock Equations
%     P*v = Rgas*T
%	  rho1*u1 = rho2*u2
%	  P1+rho1*u1^2 = P2+rho2*u2^2
%	  h1+(1/2)*u1^2 = h2+(1/2)*u2^2
% by minimizing the following pair of equations (f1, f2)
%     f1 = (P2/P1-1)+(u1^2/(P1*v1))*(v2/v1-1)
%	  f2 = (h2-h1)/((1/2)*u1^2)+(v2^2/v1^2-1)
% use Newton-Raphson Method
% Assign ideal to first guess 
    P2guess = P2ideal;
    T2guess = T2ideal; 
    v2 = v2ideal;
% NR iteration parameters
    maxTimeNR = 10; % seconds
    maxIterNR = 100; % iterations
    tolNR = 1e-10; % convergence tolerance
    stopperNR = false;
% NR loop
    iterNR = 0; 
    startTimeNR = tic;
    while stopperNR == false;
        % inc the counter
            iterNR = iterNR+1;
        % compute cp and h for gas in region 2 (see subfunction below)
            if strcmp(solnMethod,'EE')
                % compute from polynomial fits
                    [Cp2, h2] = getCph_T(T2guess, XtestGas, gasMolWt, Ru, coefVec, speciesVec); % J/kg*K, J/kg            else % 'FE' or 'FF'
            else % 'FE' or 'FF'
                % specific heat is unchanged 
                    Cp2 = Cp1; % J/kg*K
                % enthalpy change is due only to translation and rotation
                    h2 = hForm + e1V + Cp1TR*(T2guess - Tref); % J/K
            end
        % compute f1
            f1 = (P2guess/P1-1)+(u1^2/(P1*v1))*(v2/v1-1);
        % compute f2
            f2 = (h2-h1)/((1/2)*u1^2)+((v2/v1)^2-1);
        % compute derivatives for jacobian
        % note that dh/dP = 0 (ideal gas)
        % note that dT/dP = 0 (independent variables)
            df1dP2 = (1/P1)+(u1^2/(P1*v1^2))*(-v2/P2guess); % df1/dP2
            df1dT2 = (1/P1)*(0)+(u1^2/(P1*v1^2))*(v2/T2guess);
            df2dP2 = (2/u1^2)*(0)+(2*v2/v1^2)*(-v2/P2guess);
            df2dT2 = (2/u1^2)*(Cp2)+(2*v2/v1^2)*(v2/T2guess);
        % form Jacobian matrix [df1/dP2 df1/dT2; df2/dP2 df2/dT2]
            J = [df1dP2 df1dT2; ...
                 df2dP2 df2dT2];
        % matrix multiplication to solve for new T2 and P2 guesses
        % using matlab backslash operator rather than inverse saves time!
        % general method: [P2new T2new]' = [P2old T2old]' - inv(J)*[f1old f2old]'
            PTold = [P2guess T2guess]';
            fvec = [f1 f2]';
            PTnew = PTold - J\fvec;
            P2new = PTnew(1);
            T2new = PTnew(2);
        % check results - the T2 and P2 we guess should give the same T2,P2
            checkP = abs((P2guess-P2new)/P2new); % percent change from new
            checkT = abs((T2guess-T2new)/T2new);
        % decide whether to keep going or quit
            if (checkP < tolNR) && (checkT < tolNR)
                stopperNR = true; % exit the loop
            elseif toc(startTimeNR) > maxTimeNR;
                error('Timed out on NR loop')
            elseif iterNR > maxIterNR
                error('Itered out on NR loop')
            else
                % reassign the variables for the next iteration
                P2guess = P2new;
                T2guess = T2new;
                % recompute specific volume
                v2 = Rgas*T2guess/P2guess; % m3/kg
                % prevent the algorithm from going the wrong way
                while P2guess < 0
                    P2guess = mean([P2ideal P2guess]); % average back to ok range 
                end
                while T2guess < 0
                    T2guess = mean([T2ideal T2guess]);
                end
            end
    end
% save results
    P2 = P2guess;
    T2 = T2guess;
% compute gas velocity in region 2
    u2 = u1*(v2/v1); % m/s
% compute gamma in region 2
    Cv2 = Cp2 - Rgas; % J/kg*K
    gamma2 = Cp2/Cv2;
% compute density in region 2
    rho2 = 1/v2; % kg/m3
% compute speed of sound in region 2
    Vsound2 = sqrt(gamma2*P2/rho2); % m/s
% check to make sure shock conservations equations are solved
    tolCheck = 1e-6;
    eq1 = ((rho1*u1)-(rho2*u2))/(rho1*u1); % (kg/m3)*(m/s)
    eq2 = ((P1+rho1*u1^2)-(P2+rho2*u2^2))/(P1+rho1*u1^2); % Pa + (kg/m3)*(m2/s2)*(1 Pa / 1 kg/ms2)
    eq3 = ((h1+(1/2)*u1^2)-(h2+(1/2)*u2^2))/(h1+(1/2)*u1^2); % J/kg*(1 kgm2/s2 / 1 J) + m2/s2
    if abs(eq1)>tolCheck
        error('Continuity not holding over incident shock')
    elseif abs(eq2)>tolCheck
        error('Momentum not holding over incident shock')
    elseif abs(eq3)>tolCheck
        error('Energy not holding over incident shock')
    end

%% SOLVE REFLECTED SHOCK CONDITIONS

% where are we?
    %fprintf('Solving for Reflected Shock Conditions\n')
% Before Newton-Rhapson guess using ideal relations:
%     alpha=(gamma+1)/(gamma-1)
%     P5/P2=(alpha+2-P1/P2)/(1+alpha*P1/P2)
%     T5/T1=(P5/P2)*(alpha+P5/P2)/(1+alpha*P5/P2)
%	  gamma=cp1/(cp1-1) --> assumes driven gas is mostly Ar with Cp ~ 2.5    
% compute alpha
    alpha = (gamma2+1)/(gamma2-1);
% compute ideal pressure in region 2
    P5ideal = P2*(alpha+2-P1/P2)/(1+alpha*P1/P2); % Pa
% compute ideal temperature in region 2
    T5ideal = T2*(P5ideal/P2)*(alpha+P5ideal/P2)/(1+alpha*P5ideal/P2); % K
% specific volume
    v5ideal = Rgas*T5ideal/P5ideal; % m3/kg
% Solve the Reflected Shock Equations
%     P*v = Rgas*T
%     rho2*u2dash=rho5*u5
%     P2+rho2*u2dash^2=P5+rho5*u5^2
%     h2+(1/2)*u2dash^2=h5+(1/2)*u5^2
%     u2dash=u5+u1-u2
%     u2dash=(u1-u2)*(rho5/(rho5-rho2))
%     U5=0 (gas relative to lab frame behind reflected shock is stagnant)
% by solving the following pair of equations
%     f3=(P5/P2-1)+(u1-u2)^2/(P2*(v5-v2))
%     f4=(h5-h2)/((1/2)*(u1-u2)^2)+(v5+v2)/(v5-v2)
% use Newton-Raphson Method
% Assign ideal to first guess
    P5guess = P5ideal;
    T5guess = T5ideal; 
    v5 = v5ideal;
% NR iteration parameters
    maxTimeNR = 10; % seconds
    maxIterNR = 100; % iterations
    tolNR = 1e-10; % convergence tolerance
    stopperNR = false;
% NR loop
    iterNR = 0; 
    startTimeNR = tic;
    while stopperNR == false;
        % inc the counter
            iterNR = iterNR+1;
        % compute cp and h for gas in region 5 (see subfunction below)
            if strcmp(solnMethod,'EE') || strcmp(solnMethod,'FE')
                % compute from polynomial fits
                    [Cp5, h5] = getCph_T(T5guess, XtestGas, gasMolWt, Ru, coefVec, speciesVec); % J/kg*K, J/kg
            else % 'FF'
                % specific heat is unchanged 
                    Cp5 = Cp2; % J/kg*K
                % enthalpy change is due only to translation and rotation
                    h5 = hForm + e1V + Cp1TR*(T5guess - Tref); % J/K
            end
        % compute f3
            f3 = (P5guess/P2-1)+(u1-u2)^2/(P2*(v5-v2));
        % compute f4
            f4 = (h5-h2)/((1/2)*(u1-u2)^2)+(v5+v2)/(v5-v2);
        % compute derivatives for jacobian
        % note that dh/dP = 0 (ideal gas)
        % note that dT/dP = 0 (independent variables)
            df3dP5 = (1/P2)+(-(u1-u2)^2/(P2*(v5-v2)^2))*(-v5/P5guess); % df3/dP5
            df3dT5 = (1/P2)*(0)+(-(u1-u2)^2/(P2*(v5-v2)^2))*(v5/T5guess);
            df4dP5 = (1/(1/2*(u1-u2)^2))*(0)+(-2*v2/(v5-v2)^2)*(-v5/P5guess);
            df4dT5 = (1/(1/2*(u1-u2)^2))*(Cp5)+(-2*v2/(v5-v2)^2)*(v5/T5guess);
        % form Jacobian matrix [df1/dP2 df1/dT2; df2/dP2 df2/dT2]
            J = [df3dP5 df3dT5; ...
                 df4dP5 df4dT5];
        % matrix multiplication to solve for new T5 and P5 guesses
        % using matlab backslash operator rather than inverse saves time!
        % general method: [P5new T5new]' = [P5old T5old]' - inv(J)*[f3old f4old]'
            PTold = [P5guess T5guess]';
            fvec = [f3 f4]';
            PTnew = PTold - J\fvec;
            P5new = PTnew(1);
            T5new = PTnew(2);
        % check results - the T2 and P2 we guess should give the same T2,P2
            checkP = abs((P5guess-P5new)/P5new); % percent change from new
            checkT = abs((T5guess-T5new)/T5new);
        % decide whether to keep going or quit
            if (checkP < tolNR) && (checkT < tolNR)
                stopperNR = true; % exit the loop
            elseif toc(startTimeNR) > maxTimeNR;
                error('Timed out on NR loop')
            elseif iterNR > maxIterNR
                error('Itered out on NR loop')
            else
                % reassign the variables for the next iteration
                P5guess = P5new;
                T5guess = T5new;
                % recompute specific volume
                v5 = Rgas*T5guess/P5guess; % m3/kg
                % prevent the algorithm from going the wrong way
                while P5guess < 0
                    P5guess = mean([P5ideal P5guess]); % average back to ok range 
                end
                while T5guess < 0
                    T5guess = mean([T5ideal T5guess]);
                end
            end
    end
% save results
    P5 = P5guess;
    T5 = T5guess;
% compute gas velocity in region 5
% given the assumption that the gas is stagnated after the reflected shock
% in region 5, u5 (the velocity of the test gas in region 5 relative to the 
% reflected shock) is equal and opposite to the velocity of the reflected 
% shock (Urs).
    u5 = v5*(u1-u2)/(v2-v5); % m/s
    Urs = u5; % m/s 
% compute gamma in region 5
    Cv5 = Cp5 - Rgas; % J/kg*K
    gamma5 = Cp5/Cv5;
% compute density in region 5
    rho5 = 1/v5; % kg/m3
% compute speed of sound in region 5
    Vsound5 = sqrt(gamma5*P5/rho5); % m/s
% mach number of reflected shock
    u2dash = u5+u1-u2; % m/s
    Mrs = u2dash/Vsound2; % no units
% check answers
    if (T5<=T2) || (P5<=P2) || (Mrs>=Mis)
        error('Impossible reflected shock conditions')
    end
% check to make sure shock conservations equations are solved
    tolCheck = 1e-6;
    eq1 = ((rho2*u2dash)-(rho5*u5))/(rho2*u2dash); % (kg/m3)*(m/s)
    eq2 = ((P2+rho2*u2dash^2)-(P5+rho5*u5^2))/(P2+rho2*u2dash^2); % Pa + (kg/m3)*(m2/s2)*(1 Pa / 1 kg/ms2)
    eq3 = ((h2+(1/2)*u2dash^2)-(h5+(1/2)*u5^2))/(h2+(1/2)*u2dash^2); % J/kg*(1 kgm2/s2 / 1 J) + m2/s2
    if abs(eq1)>tolCheck
        error('Continuity not holding over reflected shock')
    elseif abs(eq2)>tolCheck
        error('Momentum not holding over reflected shock')
    elseif abs(eq3)>tolCheck
        error('Energy not holding over reflected shock')
    end
    
%% COMPUTE STOICHIOMETRY
    
% where are we?
    %fprintf('Solving for Stoichiometry\n')
% function to compute stoichiometry (see below)
    phi = getStoich(XtestGas, nameEle, numEle);
    
%% ASSEMBLE OUTPUT VECTORS

% where are we?
    %fprintf('Assembling Output Vectors\n')
% output vectors
    Tout = [T2 T5]; % K
    Pout = [P2 P5]; % Pa
    rhoOut = [rho1 rho2 rho5]; % kg/m3
    gammaOut = [gamma1 gamma2 gamma5]; 
    VsoundOut = [Vsound1 Vsound2 Vsound5]; % m/s
    Mout = [Mis Mrs]; % m/s
    
    


    
function specMolWts = getSpecMolecWts(nameEle, numEle, eleMolWts)
% This function computes the species molecular weights

% Find molecular weight of each species (g/mol)
    % rows of nameEle represent species
    % columns of nameEle are the elements composing that species
    specMolWts = zeros(1,size(nameEle,1)); % row vector; each entry is mol wt of the species
    gotIt = false; % used to determine if an element was found
    for i=1:size(nameEle,1) % outer loop over all rows of nameEle
        for j=1:size(nameEle,2) % inner loop over all columns of nameEle 
            tempEl = nameEle{i,j};
            if isempty(tempEl)
                continue; % if nothing there, skip it. 
            else
                for k=1:size(eleMolWts,1) %loop over rows of eleMolWts
                    if strcmp(eleMolWts{k,1},tempEl)
                        % add to the total weight; 
                        % the amount to add is the molecular weight of the 
                        % element times the number of that element present 
                        % in the species
                        specMolWts(i) = specMolWts(i)+eleMolWts{k,2}*numEle(i,j);
                        gotIt = true; % we found the element
                        break; % stop looking, we got it;
                    end
                end
                if gotIt == false
                    error('Atomic weight of %1s not found',tempEl);
                end
            end
            gotIt = false; 
        end
    end
    
function drot = getSpecDrot(speciesVec, nameEle, numEle)
% This function computes the rotational degrees of freedom of the species
    
% find rotational degrees of freedom of each species
    % list species which are linear but have more than 2 atoms
    linearSpecies = {'C3','N3','C3O2','HCN','CH2','CO2','N2O',...
        'CNC','CCN','NCO','NCN','CNN','CCO','HCC','CCCC',...
        'C2H2','C2N2','C3N2','C3H2','BECL2','F2XE','OCS','HCCCL'};
    % determine rotational degrees of freedom
    drot = zeros(1,size(nameEle,1)); % row vector
    for i=1:size(nameEle,1) % loop over all rows of nameEle
        atomTot = sum(numEle(i,:)); % add up total number of atoms in species i
        if atomTot == 1 
            drot(i) = 0; % monatomic species
        elseif atomTot == 2 
            drot(i) = 2; % linear molecule
        else % if atomTot > 2
            if sum(strcmp(speciesVec{i},linearSpecies))==1; % if the species is listed in the linear list
                drot(i) = 2; % linear molecule
            else
                drot(i) = 3; % nonlinear molecule
            end
        end

    end

    
function [Cp, h] = getCph_T(T, X, gasMolWt, Ru, coefVec, speciesVec)
% This finds the specific heat and enthalpy at the specified temperature

% loop over all species to compute for each species first
    HoverRT = zeros(1,size(coefVec,3));
    CpoverR = zeros(1,size(coefVec,3));
    for i=1:size(coefVec,3) % loop over 3rd dimension of coefVec 
        % grab temperature limits
            Tlow = coefVec(3,1,i);
            Tcom = coefVec(3,2,i);
            Thigh = coefVec(3,3,i);
        % check temperature limits
            if T < Tlow
                %fprintf('%3s: T (%1.2f K) < Tlow (%1.0f K)\n',speciesVec{i},T,Tlow);
            elseif T > Thigh
                %fprintf('%3s: T (%1.2f K) > Thigh (%1.0f K)\n',speciesVec{i},T,Thigh);
            end
        % decide which set of data to use
            if T < Tcom % low temp set
                j = 2;
            else % high temp set for T >= Tcom
                j = 1;
            end
        % grab temperature coefficients
            a1 = coefVec(j,1,i);
            a2 = coefVec(j,2,i);
            a3 = coefVec(j,3,i);
            a4 = coefVec(j,4,i);
            a5 = coefVec(j,5,i);
            a6 = coefVec(j,6,i);
            a7 = coefVec(j,7,i);
        % compute Cp/R
            CpoverR(i) = a1 + a2*T + a3*(T^2) + a4*(T^3) + a5*(T^4);
        % compute h/RT
            HoverRT(i) = a1 + a2*T/2 + a3*(T^2)/3 + a4*(T^3)/4 + a5*(T^4)/5 + a6/T;
    end
% now compute mixture properties
    Cp = (X*CpoverR')*Ru*(1000/gasMolWt); % (no units)*(J/mol*K)*(1000 g/kg)*(mol/g) = J/kg*K
    h = (X*HoverRT')*T*Ru*(1000/gasMolWt); % J/kg
 
    
function phi = getStoich(X, nameEle, numEle)
% This function computes the stoichiometry of the shock

% counting vectors
    Ccount = 0;
    Hcount = 0;
    Ocount = 0;
    Scount = 0; 
% loop over rows of nameEle (so loop over species)
    for i=1:size(nameEle,1)
        % loop over columns of nameEle (so loop over atoms composing species)
        for j=1:size(nameEle,2)
            % grab the name of the element (ie, 'H', 'C', 'He')
                tempEl = nameEle{i,j};
            % decide what to do
                if isempty(tempEl)
                    continue; % if nothing there, skip it. 
                elseif strcmp(tempEl,'C')
                    % mole fraction of species * number of atoms in species
                    Ccount = Ccount + X(i)*numEle(i,j);
                elseif strcmp(tempEl,'H')
                    Hcount = Hcount + X(i)*numEle(i,j);
                elseif strcmp(tempEl,'O')
                    Ocount = Ocount + X(i)*numEle(i,j);
                elseif strcmp(tempEl,'S')
                    Scount = Scount + X(i)*numEle(i,j);
                else
                    continue % not one of the atoms we care about
                end
        end
    end
% compute how many O atoms we need for stoichiometric burning
% C --> CO2, H --> H2O, S --> SO2
    Oneed = Ccount*2 + Hcount*(1/2) + Scount*2;
% compute phi
    phi = Oneed/Ocount; % if rich, phi>1, if lean, phi<1      
    
    
function [coefVec, nameEle, numEle] = readSandia(speciesVec, nasaFile)
% This function reads the thermodynamic data file

% See the accompanying document explaining the SANDIA.DAT file
% there are two fits for data: a high and a low temperature fit
% Each fit has 7 coefficients; we need to reassemble these. 
% Note: we are not reading in the date, and also there is a value in
% the 5th position in line 4 which is the enthalpy of formation divided
% by the gas constant R; this can be computed using the equations from
% the coefficients already read into the file.
% for more information, see the following resources:
%    http://garfield.chem.elte.hu/burcat/THERM.DAT
%    http://www.galcit.caltech.edu/EDL/public/formats/chemkin.html
%    http://www.galcit.caltech.edu/EDL/public/formats/nasaold.html
%    http://ntrs.nasa.gov/archive/nasa/casi.ntrs.nasa.gov/19780009781_1978009781.pdf

% make a checklist vector to ensure we got everything
    checkList = zeros(1,length(speciesVec));
% open the file
    fileID = fopen(nasaFile);
% read the first line
    firstLine = fgets(fileID);
% is this the right starter file?
    tempIndex = strfind(firstLine,'THERMO'); % index where 'THERMO' begins
    if isempty(tempIndex) % if 'THERMO' does not occur in first line
        error('%1s is empty',nasaFile)
    end
% read the next line
    secondLine = fgets(fileID);
% find the general temperature limits
    genTlow = str2num(secondLine(1:10)); % general low temp limit (K)
    genTcom = str2num(secondLine(11:20)); % general common temp (K)
    genThigh = str2num(secondLine(21:30)); % general high temp limit (K)
% read the file
    closeNow = false;
    while (closeNow ~= true)
        % reset variables
            elements = {''};
            numAtoms = [];
            Tlow = 0;
            Tcom = 0;
            Thigh = 0;
            TlowWarn = false;
            TcomWarn = false;
            ThighWarn = false;
        % get the first line of the species
            L1 = fgets(fileID);
        % check for ending
            endText = upper(strtrim(L1(1:3)));
            if strcmp(endText,'END')
                closeNow = true;
        % have we found everything?
            elseif sum(checkList) == length(checkList); % all 1's so this sums to the length if full! 
                closeNow = true; 
        % if not at the end, keep reading the file
            else
                % grab the species name
                    speciesName = upper(strtrim(L1(1:18))); % remove white space, make UPPERCASE
                % check the name; sometimes there is other text beside it
                    refIndex = strfind(speciesName,'REF'); % find where 'REF' occurs
                    if ~isempty(refIndex) % if 'REF' occured, then... 
                        speciesName = strtrim(speciesName(1:refIndex-1)); % take a substring
                    end
                % grab the elements and number of atoms of each
                    elpos = [25, 30, 35, 40, 74];
                    for i=1:5
                        % generate indicies
                            tmpIdx1 = elpos(i)+0; %element char 1
                            tmpIdx2 = elpos(i)+1; %element char 2
                            tmpIdx3 = elpos(i)+2; %number of elem digit 1
                            tmpIdx4 = elpos(i)+4; %number of elem digit 3
                        % grab the element name and the number of atoms
                            tmpEl = upper(strtrim(L1(tmpIdx1:tmpIdx2))); % UPPERCASE
                            tmpNum = str2num(L1(tmpIdx3:tmpIdx4));
                        % make sure the element is there 
                        % not all compounds have 5 elements
                        % if it's there, store it.  if not, throw it away. 
                        check = strcmp(tmpEl,''); % true if nothing, false if something
                        % some files have a zero in column 74
                        % this throws off the algorithm so check for it 
                        check2 = strcmp(tmpEl,'0'); 
                        if (check == false) && (check2 == false)
                            elements{1,i} = tmpEl;
                            numAtoms(1,i) = tmpNum;
                        end
                    end
                % grab the phase
                    phase = upper(L1(45)); % UPPERCASE letter
                % grab the low, middle, and high temperatures (K)
                    Tlow = str2num(L1(46:55));
                    Thigh = str2num(L1(56:65));
                    Tcom = str2num(L1(66:73));
                % check the inputs
                    if sum(size(Tlow))>2 || sum(size(Tcom))>2 || sum(size(Thigh))>2
                        error('SANDIA.DAT file format issue for %7s; Check for tabs instead of spaces.',speciesName)
                    end
                % if we did not see Tlow etc, use the general value
                    if isempty(Tlow); 
                        Tlow = genTlow;
                        TlowWarn = true; % flag that temperature was replaced
                    end
                    if isempty(Tcom);
                        Tcom = genTcom;
                        TcomWarn = true; % flag that temperature was replaced
                    end
                    if isempty(Thigh);
                        Thigh = genThigh;
                        ThighWarn = true; % flag that temperature was replaced
                    end
                % get line two
                    L2 = fgets(fileID);
                % grab the temperature coefficients here
                    Ha1 = str2num(L2(1:15)); % high temperature coefficient a1
                    Ha2 = str2num(L2(16:30));
                    Ha3 = str2num(L2(31:45));
                    Ha4 = str2num(L2(46:60));
                    Ha5 = str2num(L2(61:75));
                % get line three
                    L3 = fgets(fileID);
                % grab the temperature coefficients here
                    Ha6 = str2num(L3(1:15));
                    Ha7 = str2num(L3(16:30));
                    La1 = str2num(L3(31:45)); % low temperature coefficient a1
                    La2 = str2num(L3(46:60));
                    La3 = str2num(L3(61:75));            
                % get line four
                    L4 = fgets(fileID);
                % grab the temperature coefficients here
                    La4 = str2num(L4(1:15));
                    La5 = str2num(L4(16:30));
                    La6 = str2num(L4(31:45));
                    La7 = str2num(L4(46:60));              
                % if we want this species, store it, otherwise disregard
                    for i=1:length(speciesVec);
                        % if the same name as one of our wanted species, 
                        % and if it is a gas
                        % and if we have not grabbed it yet
                        if strcmp(speciesName,speciesVec{i}) && strcmp(phase,'G') && checkList(i)==0; 
                            % coefficeint and temperature matrix
                            coefVec(1,:,i) = [Ha1 Ha2 Ha3 Ha4 Ha5 Ha6 Ha7];
                            coefVec(2,:,i) = [La1 La2 La3 La4 La5 La6 La7];
                            coefVec(3,:,i) = [Tlow Tcom Thigh 0 0 0 0];
                            % element name and number vectors
                            for j=1:length(elements)
                                nameEle(i,j) = elements(1,j);  % use parentheses not curly brackets
                                numEle(i,j) = numAtoms(1,j);
                            end
                            % display warnings about changing temperatures
                            if TlowWarn
                                fprintf('WARNING: Defaulting Tlow to %4.2f K for %7s\n',genTlow,speciesName);
                            end
                            if TcomWarn
                                fprintf('WARNING: Defaulting Tcom to %4.2f K for %7s\n',genTcom,speciesName);
                            end
                            if ThighWarn
                                fprintf('WARNING: Defaulting Thigh to %4.2f K for %7s\n',genThigh,speciesName);
                            end
                            % check off that we got it
                            checkList(i) = 1; 
                        end
                    end
            end
    end
% close the file
    fclose(fileID);
% did we get everything?
    if sum(checkList) ~= length(checkList)
        for i=1:length(checkList)
            if checkList(i) == 0
                error('Species not in %5s: %5s',nasaFile,speciesVec{i})
            end
        end
    end
    