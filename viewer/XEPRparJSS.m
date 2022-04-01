function res = XEPRparJSS(par)

% DSC type of files
res.x=[];
res.xlabel = '';
res.y=[];
res.ylabel = '';
res.step = 0;
res.z=[];
res.zlabel = '';

% 1-st dimention ('x')

if isfield(par, 'JSS')
    switch str2num(par.JSS)
    case 2
        % cw files from Bill's Q band
        if isfield(par,'RES')
            dim = str2num(par.RES);
        else
            dim = 1024;  
        end
        center = str2num(par.HCF);
        width = str2num(par.HSW);
        res.step = width/(dim-1);
        res.x = [center-width/2:res.step:center+width/2].';
        res.xlabel = 'Magnetic field, Gs';
    case 32
        % ESP 380
        dim = str2num(par.XPLS);
        val = str2num(par.XPD9);
        res.step = val(6)*8;
        res.x = res.step*[0:dim-1].';
        res.xlabel = 'Time, ns';
    otherwise
        if isfield(par,'RES') dim = str2num(par.RES);
        elseif isfield(par,'FRES') dim = str2num(par.FRES);
        elseif isfield(par,'XPLS') dim = str2num(par.XPLS);
        else dim = 1024;
        end
        if isfield(par,'GST')
            cf = str2num(par.GST);
            wd = str2num(par.GSI);
            cf = cf + wd/2;
        else
            cf = str2num(par.HCF);
            if isfield(par, 'HSW'),
                wd = str2num(par.HSW);
            else, wd = abs(cf - str2num(par.GST))*2;
            end   
        end            
        step = wd/(dim-1);
        sf = cf - wd/2;
        res.x = sf + step * [0:dim-1].';
        if isfield(par,'JEX') res.xlabel = par.JEX;
        else res.xlabel = '?';            
        end
        if isfield(par,'JUN') res.xlabel = [res.xlabel, ', ', par.JUN];
        else res.xlabel = [res.xlabel, ', ?'];
        end
    end    
end
if isfield(par, 'MF')
    res.freq1 = str2num(par.MF)*1e9;
end

function varargout = stripunit(str)
varargout{1} = str((str>='0' & str<='9') | str=='.' | str=='E');
