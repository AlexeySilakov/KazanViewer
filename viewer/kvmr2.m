function [Freqs, Spec] = kvmr2(Sy,Ex,Op)
Parameters ={'A', 'Apa', 'Q', 'Qpa', 'Nucs', 'gn', 'I', 'nNucs'};
% lets find number of simmulaitons (parameters can be common, it means size(XX, 1)==1)
Sy.Apa = safeget(Sy, 'Apa', [0, 0, 0]);% degree
Sy.Qpa = safeget(Sy, 'Qpa', [0, 0, 0]);% degree 
Sy.nNucs = safeget(Sy, 'nNucs', 1); % number of equivalent nuclei
sim_size = 1;
isfield_num = [];
for pk = 1:length(Parameters)
    if isfield(Sy, Parameters{pk})        
        tPar = getfield(Sy, Parameters{pk});
        isfield_num(end+1) = pk;        
        if size(tPar, 1)>1
            sim_size = size(tPar, 1);
        end
    end
end
if sim_size>1
    if size(Sy.I, 2) > 1,
        Sy.I = Sy.I.';
    end
    if size(Sy.gn, 2) > 1,
        Sy.gn = Sy.gn.';
    end
    for kk = 1:length(isfield_num)
        pk = isfield_num(kk); % kind of optimisation
        tPar = getfield(Sy, Parameters{pk});
        if size(tPar, 1)==1
            Sy = setfield(Sy, Parameters{pk}, tPar(ones(sim_size, 1), :));
        elseif size(tPar, 1)~=sim_size
            error({['Parameter Sys.', Parameters{pk}, ' has wrong number of rows.']; ['Correct number is ', num2str(sim_size)]});
        end
    end
end

% Sy.Apa = safeget(Sy, 'Apa', 0);
% Sy.Qpa = safeget(Sy, 'Qpa', 0);
Sys(1:3) = Sy.g;% g-values
Sys(4:6) = 0; % angles
Sys(7) = Ex.ExciteWidth; % excitation linewidht ???? [MHz]->[cm-4]

Sys(8) = sim_size;               % number of nuclei (10 maximum)
if Sys(8) >10, error('too many nuclei (10 maximum)'); end
for ii = 1: Sys(8)
    nn = (ii-1)*14;
Sys([9:11]+nn) = Sy.A(ii, :)/2.99792458;  % HF coupling [MHz]->[cm-4]
Sys([12:14]+nn) = Sy.Apa(ii, :);             % HF tensor rotation [degree]
Sys([15:17]+nn) = Sy.Q(ii, :)/2.99792458;    % quadrupole coupling [MHz]->[cm-4]
Sys([18:20]+nn) = Sy.Qpa(ii, :);             % Q tensor rotation [degree]
Sys([21]+nn) = Sy.gn(ii);           % gn - nuclear g-value
Sys([22]+nn) = 2*Sy.I(ii)+1;                % 2*I+1
end

Exp(1) = Ex.nPoints;               % nPoints
Exp(2) = Ex.mwFreq*1e3;             % MW freq. [GHz]->[MHz]
Exp(3) = Ex.Field*10;              % Magnetic field [mT]->[G]
Exp(4) = safeget(Ex, 'tau', 0)*1e-3;                 % tau [ns]->[us]
Exp(5) = safeget(Ex, 'MaxFreq', 10);                % maximum freq. [MHz]
Exp(6) = 1;                 % xtype (gaussian)

Opt(1) = Op.nKnots;                % nKnots
Opt(2) = 0;                 % idebug

Spec = MrWinML(Sys.', Exp.', Opt.');
Spec(1) = 0;
Freqs  = linspace(0, Exp(5), Ex.nPoints);

if safeget(Sy, 'lw', 0)
Spec = convspec(Spec,1,safeget(Sy, 'lw', 0)/Exp(5)*Ex.nPoints);
end