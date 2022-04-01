% tweak pepper for 3xCu centers
function [x, Spec] = pepperoni(Sys,Exp,Opt)
if ~isempty(safeget(Sys, 'I', []))
    nSpins = length(Sys.S);
    AngFact = 1; % assume the angles are in radians
    %     [Spin,gn] = nucdata(Sys.Nucs);
    
%     if strcmp(HF_AngUn, 'rad')
%         AngFact = 1;
%     else
%         AngFact = pi/180;
%     end
    
    %     if isempty(Nucs)
    %         I = safeget(Sys, 'I', []);
    %         gn = safeget(Sys, 'gn', []); % not really necessary
    %     else
    %         [I, gn] = nucdata(Nucs); %%% to do ... make independent from easyspin
    %     end
    
    %     const = % convert MHZ to gvalues using mwFreq
    
    HF_A = safeget(Sys, 'A', [0, 0, 0]);
    if size(Sys.A, 1)~=nSpins, error('for this program size(Sys.A, 1) should be equal to size(Sys.S)'); end
    
    if size(HF_A, 1)~=length(Sys.I),
        error(['N. sets of A = ', num2str(size(HF_A, 1)), ' while length of Sys.I or Sys.Nucs = ' num2str(length(Sys.I))]);
    end
    
    HF_Apa = safeget(Sys, 'Apa', HF_A*0);
    
    for ii = 1:length(Sys.I)
        
        %%%%%%%%%%%%% not sure that this is completely right
        RM = erot(HF_Apa(ii, :)*AngFact);
        hfc = diag(HF_A(ii, :))*RM;
        hfc = sqrt(hfc(1, :).^2+hfc(2, :).^2+ hfc(3, :).^2);
        
        g_HF = Sys.g(ii, :) - Sys.g(ii, :)*Exp.mwFreq./(Exp.mwFreq+hfc*1e-3);
        
        Bs{ii} = local_splitBs(Sys.g(ii, :), Sys.I(ii), g_HF);
        
    end
    Bizz = Bs;
    nSubs = size(Bizz{1}, 1);
    
    % make pepper forget about nuclear spins
    tSys = Sys;
    tSys = rmfield(tSys, 'I');
    tSys = rmfield(tSys, 'A');
    if isfield(Sys, 'Apa'), tSys = rmfield(tSys, 'Apa'); end
    
    tExp = Exp;
        tExp.Range = Exp.Range/10; % Exp range is in gauss on the imput
    tOpt = Opt;
    Spec = zeros(1, Exp.nPoints);
    for jj = 1:nSubs
        for ii = 1:nSpins
            tSys.g(ii, :) = Bizz{ii}(jj, :);
        end
        [x, tSpec] = pepper(tSys,tExp,tOpt);
        Spec = Spec+tSpec;
    end
    
else
    Exp.Range = Exp.Range/10;
    [x, Spec] = pepper(Sys,Exp,Opt);
end
x = x*10; % Gauss is expected on the output

%%% calculate splitting
function Bs = local_splitBs(Bs, I, HF_F)

tBs = [];
for ii = 1:size(Bs, 1)
    for proj = (-I):I
        tBs(end+1, 1:3) = Bs(ii, 1:3) + proj*HF_F;
    end
end
Bs = tBs;