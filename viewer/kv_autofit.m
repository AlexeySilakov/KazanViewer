% Routines for working with formulas 
% [tokens]                 = kv_autofit(formulae)
% [funstr, funpar, funval] = kv_autofit(strFunc, cellVar)
% [modified_script, simy]  = kv_autofit(script, x, y, rangeidx {, MaxIter} )
%                                   {...} = optional
% strFunc = 'f =n*x + b*a';
% cellVar = {'n = 1'; 'b = 2'; 'a = 5'};

% --------------------------------------------------------------------
function varargout = kv_autofit(varargin) %scr_str, x, y, rangeidx)
funlist = {'fit_deer_gaus'; 'kv_expdata';'polyfit';'lshape'; 'log10'; 'exp'; 'end';'log';'cos';'sin';'abs'; './'; '.*';'.^'; ...
        '+'; '-'; '/'; '*'; '^'; ','; ''''; '='; '('; ')'; ':'; ';'; '['; ']'};

switch nargin
    % splitter of formulae
    % --------------------------------------------------------------------
case 1,
    scr_str = varargin{1};
    [name, f] = strtok(scr_str, '='); f = f(2:end);
    idx = f=='+' | f=='-' | f=='/' | f=='*' | f=='^' | f==','| f=='''';
    for k=1:size(funlist, 1)
        idxexp = strfind(f, funlist{k});
        for l = 1:size(funlist{k}, 2)
            idx(idxexp + l - 1)=1;
        end
    end
    f(idx) = ' ';
    tokens = {};
    while 1
        [tok, f] = strtok(trim(f));
        if isempty(tok), break; end;
        % Is token an x, y axis ?
        if strcmp(tok, 'x'), continue; end;
        if strcmp(tok, 'y'), continue; end;
        % Is token an i or j ?
        if strcmp(tok, 'i') | strcmp(trim(tok), 'j'), continue; end;
        % Is token a new? /alsi
        if sum(strcmp(tokens, tok)), continue; end;
        % Is token a number ?
        if isempty(str2num(tok))
            tokens{end+1} = tok;
        end
    end
    varargout{1} = tokens;
    % --------------------------------------------------------------------
    % get arguments
    % --------------------------------------------------------------------
case 2,
    scr_str = varargin{1};
    x = varargin{2};
    string = scr_str(1:end);
    cellVar = x;
    string(end+1) = ' ';
    
    % AlSi comment: funlist must be sorted by decreasing length(funlist{i})
    % otherwise additional "variables" can occure: 'acos' -> 'a   ' (if 'cos' will checked at first)
    % now it is not so important, but in future ....
    
    tempstr = string;
    for j = 1:length(funlist)
        pos = findstr(tempstr, funlist{j});
        if isempty(pos), continue; end
        for i = 1: length(pos)
            len = length(funlist(j));
            for q = (1:len)
                tempstr(q + pos-1) = ' ';
            end
        end
    end
    shiftpos = 0;
    
    space = ' ';
    funpar = {};
    for k = 1:length(cellVar)
        
        varname = strtokstr(cellVar{k}, ' ='); % name of variable
        funpar{end+1} = varname;
        varlen = length(varname);
        
        varpos = 0;
        varpos = findstr(tempstr, [' ', varname, ' ']);
        while ~isempty(varpos)
            varstr = ['in(' num2str(k), ')']; % new variable string
            stpos = varpos(1);
            endpos = varpos(1)+varlen + 1;
            
            str1 = string(1:stpos);
            str2 = string(endpos:end);
            ss1 = tempstr(1:stpos);
            ss2 = tempstr(endpos:end);
            string = [str1, varstr, str2];
            tempstr = [ss1, space(ones(1, length(varstr))) ,ss2];
            
            eval([cellVar{k}, ';']);
            eval([varstr, '=', varname, ';']);
            
            varpos = findstr(tempstr, [' ', varname, ' ']);
        end
    end
    varargout{3} = in;
    varargout{2} = funpar;
    varargout{1} = string;
    % --------------------------------------------------------------------
    % automatic fitting
    % --------------------------------------------------------------------
otherwise
    %%%%%%%% construct fitting function %%%%%%%
    scr_str = varargin{1};
    x = varargin{2};
    y = varargin{3};
    rangeidx = varargin{4};
    if nargin ==5, 
        MaxIter = varargin{5};
    else 
        MaxIter   = 500;
    end
    var_idx = [];
    form_par = {};
    add_par  = {};
    
    for k = 1:length(scr_str)
        tempstr = scr_str{k};
        if strcmp(tempstr(1:3), 'f ='), 
            form_par{end+1} = scr_str{k};
            funstr = scr_str{k}; 
        elseif strcmp(tempstr(1:min(6,end)), 'FitRng')   %  common
            form_par{end+1} = scr_str{k};
        elseif strcmp(tempstr(1:min(7,end)), 'MaxIter')  % general fit
            eval([tempstr,';']);
        elseif strcmp(tempstr(1:min(3,end)), 'pol')      % polinomial
        else
            var_idx(end + 1) = k;     
        end
    end
    if isempty(funstr), msgbox('can not find function');  return; end
    
    % recognising variables names
    [funstr, pars, in] = kv_autofit(funstr, {scr_str{var_idx}});
    
    if isempty(rangeidx), rangeidx = [1:length(x)]; end

    %%%%%%%%%%%%% Fitting %%%%%%%%%%%%%%
    sz2 = size(y, 2);
    resstr = funstr(4:end);
    if findstr(resstr, 'exppolyfit')
        for kk = 1:sz2
            shift      = min(y(:,kk))-.0001*abs(max(y(:,kk))-min(y(:,kk)));
            pol        = polyfit(x(rangeidx), log(y(rangeidx,kk)-shift), in(1));
            simy(:,kk) = exp(polyval(pol, x)) + shift; 
        end

        for i = 1:length(pol),
            add_par{end+1} = ['pol', num2str(i), ' = ', num2str(pol(i))];
        end
    elseif findstr(resstr, 'polyfit')
        for kk = 1:sz2
            pol        = polyfit(x(rangeidx), y(rangeidx,kk), in(1));
            simy(:,kk) = polyval(pol, x); 
        end
        for i = 1:length(pol),
            add_par{end+1} = ['pol', num2str(i), ' = ', num2str(pol(i))];
        end
    else
        strChi = ['sum((', resstr,' - real(y)).^2)/1e-6'];
        f = inline(strChi, 'in', 'x', 'y');
        
        opt = [];
        opt = optimset('Display', 'off', 'MaxIter', MaxIter);
        
        for kk = 1:sz2
            in = fminsearch(f, in, opt, x(rangeidx), real(y(rangeidx, kk)));
            eval(['simy(:,kk) = ', resstr, ';']);
        end
    end
    
    % parameter string construction 
    for m = 1:length(in)
        form_par{end+1} = eval(['[pars{m},'' = '', num2str(in(', num2str(m),'))];']);
    end

    for m = 1:length(add_par)
        form_par{end+1} = add_par{m};
    end
    
    varargout{2} = simy;
    varargout{1} = form_par;
end