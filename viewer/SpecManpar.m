function res = SpecManpar(par)

% KAZAN dataviewer with plugins By Boris Epel & Alexey Silakov
% MPI of Bioinorganic Chemistry, Muelhaim an der Ruhr, 2003
% Free for non-commercial use. Use the program at your own risk. 
% The authors retains all rights. 
% Contact: epel@mpi-muelheim.mpg.de
prefix = ['n', 'u', 'm', ' ', 'k', 'M', 'G'];
koeff  = [1E-9, 1E-6, 1E-3, 1, 1E3, 1E6, 1E9];
res.title = safeget(par, 'general_name', '?');
sweepax = {};
key = 'transient';
fullfield = ['sweep_', key];
idx = 0;
triggers = str2num(safeget(par, 'streams_triggers', '1'));
while isfield(par, fullfield)
    [ax.t, str] = strtok(getfield(par, fullfield), ',');
    ax.t = trim(ax.t);
    [ax.size, str] = strtok(str(2:end), ',');
    if ax.t=='S' | ax.t=='I', ax.size = 1; else ax.size = str2num(ax.size); end
    [ax.reps, str] = strtok(str(2:end), ',');
    ax.reps = str2num(ax.reps);
    ax.var = {};
    while ~isempty(str)
        [ax.var{end+1}, str] = strtok(str(2:end), ',');
    end
    sweepax{end+1,1} = ax;
    fullfield = ['sweep_sweep', num2str(idx)];
    idx = idx +1;
end
sweepax{1}.size=sweepax{1}.size*triggers;
res.sweepax = sweepax;

axislabel = ['xyz'];
counter = 1;
for k = 1:size(sweepax, 1)
    arr = [];
    asize = sweepax{k}.size;
    if asize > 1
        if sweepax{k}.t == 'I'
            tempparam = 'trans';
            par.params_trans = '1sl step 1sl;';
        elseif sweepax{k}.t == 'T'
            tempparam = 'trans';
            par.params_trans = '1 step 1;';
        else
            tempparam = sweepax{k}.var{1};
            tempparam(findstr(tempparam, ' ')) = '_';
            tempparam(findstr(tempparam, '.')) = '_';
            tempparam(findstr(tempparam, ',')) = '_';
        end

        % check if this is a parameter
        if isfield(par, ['params_', tempparam])

            str = getfield(par, ['params_', tempparam]);
            [tk1, str1] = gettoken(str, 'to'); 
            if isempty(tk1)
                [tk1, str1] = gettoken(str, 'step'); 
                if isempty(tk1)
                    % string of the type 10ns, 20ns, 30ns;
                    [str1] = gettoken(str, ';');
                    [tk1, str1] = gettoken(str1, ','); 
                    while ~isempty(tk1)
                        [arr(end+1),unit] = kvgetvalue(tk1);
                        [tk1, str1] = gettoken(str1, ','); 
                        if isempty(tk1) & ~isempty(str1)
                            tk1 = str1; str1 = [];
                        end
                    end
                else
                    % string of the type 10ns step 6 ns
                    tk2 = trim(gettoken(str1, ';'));
                    [minval, unit] = kvgetvalue(tk1);
                    step = kvgetvalue(tk2);
                    arr = [0:asize-1]*step+minval;
                end
            else
                % string of the type 10ns to 60 ns
                tk2 = trim(gettoken(str1, ';'));
                [minval, unit] = kvgetvalue(tk1, 0);
                [maxval, maunit] = kvgetvalue(tk2, 0);
                    
                arr = [0:1/(asize-1):1]*(maxval-minval)+minval;
            end
        else
            str = getfield(par, ['aquisition_', tempparam]);
            arr = [0:asize-1];
            unit = 's';
        end
        
        % Unit normalization
        
%         switch unit
%         case 'cm^{}'
%             unit = 'cm^{-1}';
%         case 'G',
%         case 'K',
%             case 's',
%             otherwise
%                 umax = max(abs(arr));
%                 for kk = length(koeff):-1:1
%                     if umax > koeff(kk)
%                         uk = koeff(kk);
%                         unit = [prefix(kk), unit];
%                         arr = arr./uk;
%                         break;
%                     end
%                 end
%         end

        res = setfield(res, axislabel(counter), arr');
        res = setfield(res, [axislabel(counter), 'label'], ...
            [sweepax{k}.var{1}, ', ',unit]);
        counter = counter + 1;
    end
end

function [tk,rstr] = gettoken(istr, tok)

pos = strfind(istr, tok);

if isempty(pos)
    tk=[];
    rstr=istr;
else
    tk=trim(istr(1:pos-1));
    rstr=trim(istr(pos+length(tok):end));
end
