% function [outy, outax] = processing(y, inax)
% ax.filt       <'none'>/'rc'/'mvgaver'/'savgol'
% ax.filtpar1   <'auto'>/'1D'/'2D'/'1D2D'

function [outy, outax] = processing(y, inax)
if nargin < 2
    disp('Usage: [outy, outax] = processing(y, inax)');
    disp('Filter:')
    disp('ax.filt =      <''none''>/''rc''/''mvgaver''/''savgol''');
    disp('ax.filtpar1 =  <''auto''>/''1D''/''2D''/''1D2D''');
    disp('ax.tc = 1..n');
    disp('Pseudo modulation: ax.diff = <''diff''> and ax.ps_mod [units of inax.x]')
    error('Wrong syntax.');
end
if isempty(y), 
    outy=[]; outax=[]; 
    return; 
end;
filtproc = safeget(inax, 'filt', '');
diffintproc = safeget(inax, 'diff', '');
shift = safeget(inax,'s', 0);
filtx = 0; filty = 0;
switch filtproc
case 'rc',
    sizey = size(y,1);
    cc = safeget(inax, 'tc', 5);
    y = rcfilt(y, 1, cc);
case 'mvgaver',
    cc = safeget(inax, 'tc', 5);
    fpar = safeget(inax, 'filtpar1', 'auto');
    if strcmp(fpar,'auto')
        filtx = 1;
        if size(y,2) > 5, filty = 1; end
    else
        if strfind(fpar, '1D'), filtx = 1; end
        if strfind(fpar, '2D'), filty = 1; end
    end
    if filtx, y = kv_mvgavg(y, cc, 'binom'); end
    if filty, y = kv_mvgavg(y.', cc, 'binom').'; end
case 'savgol',
    cc   = safeget(inax, 'tc', 5);
    cc1  = safeget(inax, 'pol', 2);
    fpar = safeget(inax, 'filtpar1', 'auto');
    if strcmp(fpar,'auto')
        filtx = 1;
        if size(y,2) > 5, filty = 1; end
    else
        if strfind(fpar, '1D'), filtx = 1; end
        if strfind(fpar, '2D'), filty = 1; end
    end
    if filtx, y = kv_mvgavg(y, cc, 'savgol', cc1); end
    if filty, y = kv_mvgavg(y.', cc, 'savgol', cc1).'; end
otherwise
end      
switch diffintproc
case 'diff0',
    psharm = safeget(inax, 'ps_mod_harm', 1);
    for ii =1:psharm
        y = [0; diff(y)];
    end
case 'diff',
    sizey = size(y,1);
    % alsi 30.08.2009 change to the KV version of "pseumod"
%     bl = sum(y([1:4,end-3:end],:))/8;
    
    cc = safeget(inax, 'ps_mod', 5);
    psharm = safeget(inax, 'ps_mod_harm', 1);
    if cc ~= 0, 
        y = kv_pseumod(inax.x, y, cc, psharm);
    else
        for ii =1:psharm
            y = [zeros(1, size(y, 2)); diff(y, 1, 1)];
        end
    end
%     for k = 1:size(y,2)
%         y(:,k) = pseumod([1:sizey]', (y(:,k)-bl(1,k)).', cc).';
%     end
case 'integr'
    % alsi 07.01.2004  
    sizey = size(y,1);
    bl = sum(y([1:4,end-3:end],:))/8;
    sizey = size(y, 2);
    % a bit rough but enougth for viewing
    for k = 1:sizey
        y(:, k) = cumsum(y(:, k)-bl(1,k));
    end
    y = y.*(inax.x(2) - inax.x(1)); % /2
otherwise
end  
outax = inax;
outax.filt = '';
outax.diff = '';
outy = y + shift;      
