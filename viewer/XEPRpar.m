function res = XEPRpar(par)

if ~isstruct(par), error('Input must be a structure of string parameters'); end;

res = [];
% extracting some useful parameters from file
if isfield(par, 'MWFQ')
    res.freq1 = str2num(par.MWFQ);
end
if isfield(par, 'MF')
    res.freq1 = str2num(par.MF);
end
if isfield(par, 'CenterField')
    res.cf = str2num(stripunit(par.CenterField))*1E-4;
end
if isfield(par, 'IKKF')
    res.complex = strcmp(par.IKKF, 'CPLX');
elseif isfield(par, 'XQAD')
    res.complex = strcmp(par.XQAD, 'ON');
else, res.complex = 0; 
end
if isfield(par, 'TITL')
    res.title = par.TITL;
elseif isfield(par, 'JCO')
    res.title = par.JCO;
else
    res.title = '?';
end

res = getrespar321(par, res, 'x');
res = getrespar321(par, res, 'y');

function varargout = stripunit(str)
varargout{1} = str((str>='0' & str<='9') | str=='.' | str=='E');

function res = getrespar321(par, res, axletter)

if ~isstruct(par), error('Input must be a structure of string parameters'); end;

res = setfield(res, lower(axletter), []);
res = setfield(res, [lower(axletter) 'label'], '');

% defaults 
dim = 1; sf = 0; step = 1; label = '?';
unit = safeget(par, [upper(axletter) 'UNI'], '?');
unit = unit(unit~='''');

% XSophe stile 
if ~isfield(par,'JSS') & (isfield(par,'RES')|isfield(par,'RRES')) & isfield(par, 'HCF')
    par.JSS = '2';
    par.RES = safeget(par, 'RES', safeget(par,'RRES','1024')); % cheating
end

% attempt to produce valid sf, dim, step, label
if isfield(par, 'JSS')
    switch str2num(par.JSS)
    case 2
        % CW fiels sweep files 
        if ~strcmp(upper(axletter), 'Y');
            dim = str2num(safeget(par,'RES','1024'));
            if isfield(par, 'GST')& isfield(par, 'GSI')
                sf = str2num(par.GST);
                step = str2num(par.GSI)/(dim-1);
            elseif isfield(par, 'HCF')& isfield(par, 'HSW')
                center = str2num(par.HCF);
                width = str2num(par.HSW);
                step = width/(dim-1);
                sf = center-width/2;
            elseif isfield(par, 'HCF')& isfield(par, 'GST')
                center = str2num(par.HCF);
                width = 2*abs(center - str2num(par.GST));
                step = width/(dim-1);
                sf = center-width/2;
            elseif isfield(par, 'HCF')& isfield(par, 'GSI')
                par.DOS = '1';
                center = str2num(par.HCF);
                width = str2num(par.GSI);
                step = width/(dim-1);
                sf = center-width/2;                
            else, 
                width = abs(center - str2num(par.GST))*2;
                step = width/(dim-1);
                sf = center-width/2;
            end      
            label = 'Magnetic Field, G';
        end
    case 32
        % ESP 380
        if ~strcmp(upper(axletter), 'Y')
            if isfield(par, [upper(axletter) 'QNT'])
                dim = str2num(getfield(par, [upper(axletter),'PLS']));
                idx = getfield(par,[upper(axletter) 'QNT']);
                idx = idx(idx~=' ');
                switch idx
                case 'Time'
                    step = 0;
                    val = str2num(par.Psd5);
                    for i=69:75
                        if val(i) > 0
                            res.step = val(i); 
                            break;
                        end
                    end
                case 'Magn.Field'
                    if isfield(par,'HCF')
                        cf = str2num(par.HCF);
                        if isfield(par, 'HSW'),
                            wd = str2num(par.HSW);
                        else, wd = abs(cf - str2num(par.GST))*2;
                        end   
                    else
                        cf = str2num(par.GST);
                        wd = str2num(par.GSI);
                        cf = cf + wd/2;
                    end
                    step = wd/(dim-1);
                    sf = cf - wd/2;
                    label = 'Magnetic Field, G';
                case {'RF1', 'RF2'}
                    sf = sscanf(getfield(par, [idx,'StartFreq']), '%f');
                    wd = sscanf(getfield(par, [idx,'SweepWidth']), '%f');
                    step = wd/(dim-1);
                    label = [idx, ', MHz'];
                case {'1.RFSource'}
                    sf = str2num(getfield(par, 'ESF'));
                    wd = str2num(getfield(par, 'ESW'));
                    step = wd/(dim-1);
                    label = ['RF, MHz'];
                end
            else
                dim = str2num(getfield(par, [upper(axletter),'PLS']));
                if isfield(par,'JUN')
                    unit = par.JUN;
                end
                if isfield(par,'GST')
                    sf=str2num(par.GST);
                    wd = str2num(par.GSI);
                    step = wd/(dim-1);
                    label = ['?,',unit];
                else 
                    if isfield(par, 'HCF')
                        cf = str2num(par.HCF)*1E-4;
                    end        
                    dim = str2num(getfield(par, [upper(axletter),'PLS']));
                    val = str2num(par.XPD9);
                    step = val(6)*8;
                    label = 'Time, ns';
                end
            end
        end
    % 2D files         
    case {4128,4144}
        dim = str2num(getfield(par, ['SS',upper(axletter)]));        
        if strcmp(upper(axletter), 'X') & res.complex,
            dim = dim / 2;
        end
        sw = str2num(safeget(par, ['X',upper(axletter) 'WI'], num2str(dim)));
        step = sw./(dim-1);
        unit = safeget(par, ['X',upper(axletter) 'UN'], '?');
        unit = unit(unit~='''');
        label = ['Time, ', unit];
    otherwise
        if isfield(par, [upper(axletter) 'TYP'])
            % Latest Bruker format 
            if ~strcmp(getfield(par,[upper(axletter) 'TYP']), 'NODATA'), 
                nam = safeget(par, [upper(axletter) 'NAM'], '?');
                dim = str2num(getfield(par, [upper(axletter),'PTS']));
                sf = str2num(getfield(par, [upper(axletter),'MIN']));
                wd = str2num(getfield(par, [upper(axletter),'WID']));
                step = wd/(dim-1);
                label = [nam, ' ', unit];
            end
        elseif isfield(par, [upper(axletter) 'NAM'])
            % predefined experiments in ESP580
            dim = str2num(getfield(par, [upper(axletter),'PTS']));
            sf = str2num(getfield(par, [upper(axletter),'MIN']));
            wd = str2num(getfield(par, [upper(axletter),'WID']));
            step = wd/(dim-1);
            nam = safeget(par, [upper(axletter) 'NAM'], '?');
            label = [nam, ', ', unit];
        end
    end
    %   no JSS ESP580 format  
elseif isfield(par, 'JEX')
    if ~strcmp(upper(axletter), 'Y');
        dim = str2num(safeget(par,'RES','1024'));
        switch par.JEX
            case 'ENDOR'
                sf = str2num(par.ESF);
                width = str2num(par.ESW);
                step  = width/(dim-1);
                label = ['Frequency, MHz'];
        end
    end
elseif ~strcmp(getfield(par,[upper(axletter) 'TYP']), 'NODATA'),
    dim = str2num(getfield(par, [upper(axletter),'PTS']));
    sf = str2num(getfield(par, [upper(axletter),'MIN']));
    wd = str2num(getfield(par, [upper(axletter),'WID']));
    step = wd/(dim-1);
    nam = safeget(par, [upper(axletter) 'NAM'], '?');
    label = [nam, ', ', unit];
    
    if isfield(par, [upper(axletter) 'AxisQuant'])
        % advance experiment in ESP580
        idx = getfield(par,[upper(axletter) 'AxisQuant']);
        idx = idx(idx~=' ');
        label = [idx, ', ', unit];
    end
end
res = setfield(res, lower(axletter), sf + step * [0:dim-1].');
res = setfield(res, [lower(axletter) 'label'], label);

