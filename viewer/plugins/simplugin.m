function varargout = simplugin(varargin)
% SIMPLUGIN Application M-file for simplugin.fig
%    Intended to use only with 'kazan' viewer
% KAZAN dataviewer with plugins By Boris Epel & Alexey Silakov
% MPI of Bioinorganic Chemistry, Muelheim an der Ruhr, 2003
% Free for non-commercial use. Use the program at your own risk.
% The authors retains all rights.
% Contact: epel@mpi-muelheim.mpg.de

% Last Modified by GUIDE v2.5 08-Aug-2005 12:11:26
% boep 19-dec-2004, MPI

if (nargin == 0)  % LAUNCH GUI

    token = 'SIMPLUGIN_open';
    oldfig = getappdata(0, token);
    oldfig = oldfig(ishandle(oldfig));
    
    [fpath,~,e] = fileparts(which('kazan'));
    inifilename = [fpath , filesep, 'kazan.ini'];
    ini = inimanaging(inifilename);
    opc = 2;
    try
        opc = ini.KazanViewer.MultyKazan;
    catch
        opc = 2;
    end
    if ((opc == 1|| opc == 2) && ~isempty(oldfig)), 
        set(oldfig(1),'Visible','on'); varargout{1}=oldfig(1); return; 
    end
    
    if ~isempty(dir([fpath, filesep, 'simulations']))
        addpath([fpath, filesep, 'simulations']);
    end
    fig = figure; %openfig(mfilename,'new');
    setappdata(0, token, [oldfig, fig]);
    
    set(fig, 'Units', 'pixels', 'Name', 'Simplugin', 'NumberTitle', 'off');
    vv = get( 0, 'ScreenSize');
    set(fig, 'Position', [319 100 250 min(800, vv(4)-200)], 'Tag', 'figure1', 'Visible', 'off', 'MenuBar', 'none'); %%% hide until everything is done
    set(fig, 'ResizeFcn', 'simplugin(''figure1_ResizeFcn'',[],[],guidata(gcf))');
    
    % Generate a structure of handles to pass to callbacks, and store it.
    handles = guihandles(fig);
    handles.figure1 = fig;
    handles = figuresetup(handles);
    figure1_ResizeFcn(handles.figure1, [], handles)
    [pathstr, ~, ~] = fileparts(get(fig, 'FileName'));
    outstr = inimanaging([pathstr, filesep, 'kazan.ini']);
    try
        pp = outstr.Common_For_Plugins;
        handles.currdir = pp.currDir;
    catch
        handles.currdir = pwd;
    end
    try set(handles.mShsumm, 'Checked', outstr.Simplugin.showSum); end
    try set(handles.mShdiff, 'Checked', outstr.Simplugin.showDiff); end
    try set(handles.mSimAutoBaseline, 'Checked', outstr.Simplugin.autoBaseline); end
    try set(handles.mSimAgr, 'Checked', outstr.Simplugin.applSourceProc);end
    try set(handles.mSimApplyFreq, 'Checked', outstr.Simplugin.loadExpFreq); end
    guidata(fig, handles);
    set(fig,  'Visible', 'on');
    if nargout > 0,
        varargout{1} = fig;
    end
else
  if ischar(varargin{1}) % INVOKE NAMED SUBFUNCTION OR CALLBACK
    try
        %     [varargout{1:nargout}] =
        feval(varargin{:}); % FEVAL switchyard
    catch
        disp(lasterr);
    end
  end
end


    
%| ABOUT CALLBACKS:
%| GUIDE automatically appends subfunction prototypes to this file, and 
%| sets objects' callback properties to call them through the FEVAL 
%| switchyard above. This comment describes that mechanism.
%|
%| Each callback subfunction declaration has the following form:
%| <SUBFUNCTION_NAME>(H, EVENTDATA, HANDLES, VARARGIN)
%|
%| The subfunction name is composed using the object's Tag and the 
%| callback type separated by '_', e.g. 'slider2_Callback',
%| 'figure1_CloseRequestFcn', 'axis1_ButtondownFcn'.
%|
%| H is the callback object's handle (obtained using GCBO).
%|
%| EVENTDATA is empty, but reserved for future use.
%|
%| HANDLES is a structure containing handles of components in GUI using
%| tags as fieldnames, e.g. handles.figure1, handles.slider2. This
%| structure is created at GUI startup using GUIHANDLES and stored in
%| the figure's application data using GUIDATA. A copy of the structure
%| is passed to each callback.  You can store additional information in
%| this structure at GUI startup, and you can change the structure
%| during callbacks.  Call guidata(h, handles) after changing your
%| copy to replace the stored original so that subsequent callbacks see
%| the updates. Type "help guihandles" and "help guidata" for more
%| information.
%|
%| VARARGIN contains any extra arguments you have passed to the
%| callback. Specify the extra arguments by editing the callback
%| property in the inspector. By default, GUIDE sets the property to:
%| <MFILENAME>('<SUBFUNCTION_NAME>', gcbo, [], guidata(gcbo))
%| Add any extra arguments after the last argument, before the final
%| closing parenthesis.
% --------------------------------------------------------------------

function FileUpdate(handles)

handles = SetRange(handles);

if handles.applyfreq
    hand = guidata(handles.hh);
    
    % Loading of the experiment mw frequncy
    newval = safeget(hand.src.ax, 'freq1', 9.34)*1E-9;
    Str = get(handles.list,'String');
    Str = UpdatePars(Str, 'Exp.mwFreq', newval);
    set(handles.list,'String', Str);
    
    for kk=1:size(handles.sim,2)
        handles.sim{kk}.Script = UpdatePars(Str, 'Exp.mwFreq', newval);
    end
    guidata(handles.figure1, handles);
    list_Callback(handles.figure1, [], handles);
end

% --------------------------------------------------------------------
function varargout = mFileSave(h, eventdata, handles, varargin)
handles = guidata(handles.figure1);
[pathstr,name,ext] = fileparts(get(handles.figure1, 'FileName'));
handles.currdir = safeget(handles, 'currdir', pathstr);
[fname,pname] = uiputfile([handles.currdir, '*.*'],'Save Script');
handles.currdir = pname;
if fname == 0, return; end

y = real(handles.y);
x = handles.ax.x + handles.ax.dx;

p = [x,y];
save([pname, fname], 'p', '-ascii')
guidata(handles.figure1, handles);
% --------------------------------------------------------------------

function varargout = list_Callback(h, eventdata, handles, varargin)
handles = guidata(handles.figure1);
Res = get(handles.list, 'String');
val = get(handles.list, 'Value');

[name, str1] = strtok(Res{val}, '=');
str2 = [trim(str1(3:end)), ' '];

if strcmp(trim(name(3:end)), 'f')
    set(handles.ePars, 'String', str2);
    set(handles.ePars, 'UserData', 'function');
elseif str2(1)==''''
    primes = findstr(str2, '''');
    if size(primes, 2) < 2, primes(2) = size(str2, 2); end
    set(handles.ePars, 'String', str2(2:(primes(2)-1)));
    set(handles.ePars, 'UserData', 'string');
else
    dig = str2num(str2);
    set(handles.ePars, 'String', str2);
    set(handles.ePars, 'UserData', 'value');
end

% --------------------------------------------------------------------
function varargout = ePars_Callback(h, eventdata, handles, varargin)
handles = guidata(handles.figure1);
Res = get(handles.list, 'String');
val = get(handles.list, 'Value');
newval = trim(get(handles.ePars, 'String'));

[name, str1] = strtok(Res{val}, '=');
name = trim(name); opt_mode=name(1:2); name=name(3:end);
commOK = 1;
switch get(handles.ePars, 'UserData')
case 'function'
    set(handles.ePars, 'String', newval);
    %tokens = fsplitter([name ' = ' newval]);
    tokens = kv_autofit([name ' = ' newval]);
    % load script 
    for k=1:size(Res,1),
        if ~strcmp(trim(strtok(Res{k},'=')),'f')
            try,eval([Res{k}(3:end) ';']);end
        end
    end
    Res = {[opt_mode, name ' = ' newval]};
    for k = 1: size(tokens, 2)
        try, var = eval(tokens{k}); catch, var = 0; end
        Res{end+1} = [opt_mode, tokens{k},' = ', num2str(var)];
    end
    commOK = 0;
case 'string'
    set(handles.ePars, 'String', newval);
    Res{val} = [opt_mode, name ' = ''', newval,''''];
otherwise
    nums = str2num(newval);
    set(handles.ePars, 'String', trim(newval));
    Res{val} = [opt_mode, name ' = ' trim(newval)];
end
set(handles.list, 'String', Res);
list_Callback(h, eventdata, handles);
handles.sim{get(handles.popSel, 'Value')}.Script = Res;
guidata(handles.figure1, handles);
if strcmp(opt_mode(1), 'c')&commOK,
    handles = change_comm(handles, Res{val});
end
guidata(handles.figure1, handles);
if isempty(strfind(name,'Shift.'))
    if strcmp(get(handles.mCalc, 'Checked'), 'on')|(h==handles.pbReCalc), 
        nn = get(handles.popSel, 'Value');
        switch handles.sim{nn}.stype
        case 'sim'
            handles = pepcalc(get(handles.popSel, 'Value'), handles);
            handles = SetRange(handles);
        case 'baseline'
            handles = bscalc(get(handles.popSel, 'Value'), handles);
        end
        plothis(handles);
    end
else
    plothis(handles);
end

% % --------------------------------------------------------------------
% function varargout = slPars_Callback(h, eventdata, handles, varargin)
% handles = guidata(handles.figure1);
% if ~strcmp(get(handles.ePars, 'UserData'),'value'), return; end;
% shift = get(handles.slPars, 'Value');
% set(handles.slPars, 'Value', 0.5);
% num = get(handles.popup, 'Value');
% str = get(handles.popup, 'String');
% dim = str2num([str{num}]);
% 
% CurVal = str2num(get(handles.ePars, 'String'));
% Val = CurVal + dim*(shift - 0.5)*2;
% 
% set(handles.ePars, 'String', num2str(Val, 15));
% ePars_Callback(h, eventdata, handles, varargin);

% --------------------------------------------------------------------
function handles = bscalc(num, handles)
handles = guidata(handles.figure1);
hand = guidata(handles.hh);
x = hand.src.ax.x;
func = ';';
% load script 
str = get(handles.list, 'String');
for k=1:size(str,1),
    if strcmp(trim(strtok(str{k}(3:end),'=')),'f')
        func = str{k}(3:end);
    else
        try,eval([str{k}(3:end) ';']);end
    end
end
try,eval([func ';']); 
catch 
    disp(lasterr);   
end

handles.sim{num}.x = x;
handles.sim{num}.amp = abs(max(f)-min(f));
handles.sim{num}.y = f(:, ones(size(hand.src.y, 2), 1));
handles.sim{num}.MultSize = size(hand.src.y, 2);
guidata(handles.figure1, handles);

% --------------------------------------------------------------------
function mSimShow_Callback(hObject, eventdata, handles)
handles = guidata(handles.figure1);
c = {'off', 'on'};
stat = strcmp(get(hObject, 'Checked'), 'on');
set(hObject, 'Checked', [c{~stat + 1}]);
plothis(handles);

% --------------------------------------------------------------------
function varargout = pepcalc(num, handles)
handles = guidata(handles.figure1);
try
    set(handles.figure1,'name','Busy')
    % load script 
    Sys = []; Exp = []; Opt = []; Shift = [];
    %str = get(handles.list, 'String');
    str = handles.sim{num}.Script;
    for k=1:size(str,1),
        try,eval([str{k}(3:end) ';']);
        catch, disp(['Warning: Problem with string ', str{k}(3:end)]);
        end
    end
    
    hand = guidata(handles.hh);
    % Default values 
    ll = length(hand.src.ax.x);
    if ll < 5, ll=1024; end; 
    Exp.nPoints = safeget(Exp, 'nPoints', ll);
    Opt.Verbosity = safeget(Opt, 'Verbosity', 0);
    
    sim = safeget(Opt, 'Sim', 'pepper');
    
    if safeget(Exp, 'gField', 0) > .5
        Exp.Field = g2fld(Exp.gField, Exp.mwFreq*1E9)*1E3;
    end 
    %%%%%%%%%%%%%%%%%%%%%%%%%%% MULTIPLE PARAMETER FIELDS  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    MultFld1 = {'Exp', 'Exp', 'Opt', 'Exp','Sys', 'Shift', 'Exp', 'Exp', 'Sys', ...
        'Exp', 'Exp', 'Sys', 'Exp', 'Sys'};
    MultFld2 = {'Field', 'gField', 'ThetaRange', 'ExciteWidth','lwEndor', 'Scale', 'mwFreq', 'asym', 'lwEndor', ...
        'TrFreq', 'Hole', 'RelaxPar', 'Temp', 'Omega'};
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    MultSize = 1;
    % Determine multiple argument variables 
    for kk=1:size(MultFld1, 2)
        if eval(['isfield(',MultFld1{kk},',''',MultFld2{kk},''')'])
            eval([MultFld2{kk},'=',MultFld1{kk},'.',MultFld2{kk},';']);
            MultSize = eval(['max([MultSize, size(',MultFld2{kk},',1)]);']);
        end
    end
    
    % Construct correct-dimensioned variables
    for kk=1:size(MultFld1, 2)
        if eval(['exist(''',MultFld2{kk},''')'])==1
            eval(['sz = size(',MultFld2{kk},',1);']);
            if sz < MultSize
                eval([MultFld2{kk},'(',num2str(sz+1),':',num2str(MultSize),',:)=',...
                        MultFld2{kk},'(zeros(',num2str(MultSize-sz),',1)+end,:);']);
            end
        end
    end
    
    % unit conversion
    if strcmp(sim, 'pepper') | strcmp(sim, 'simpleEPR') | strcmp(sim, 'EPRsim_ext') , defltunit = 10; else defltunit = 1; end % | strcmp(sim, 'EPRpolarized')
    Opt.unit_cnv = safeget(Opt, 'unit_cnv', defltunit); 

    % save field to use later for shifting of data
    if isfield(Sys, 'Nucs') & exist('Field')
        handles.Larmor = fld2nfreq(Field/Opt.unit_cnv*1E-3, Sys.Nucs)*1E-6; % in MHz
    elseif isfield(Sys, 'gn') & exist('Field')
        handles.Larmor = Sys.gn(1)*Field/Opt.unit_cnv*1E-9*nmagn/planck; % in MHz
    elseif exist('Larmor'), handles = rmfield(handles,'Larmor'); 
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
%%%%%%%%%%%%%%%%%%%%% Simulations %%%%%%%%%%%%%%%%%%%%%%%%%     
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
    handles.is2D = ~isempty(strfind(sim, '2D'));
    for kk=1:MultSize
        Opt.SimIdx = kk;
        % Determine multiple argument variables 
        for kk1=1:size(MultFld1, 2)
            if eval(['isfield(',MultFld1{kk1},',''',MultFld2{kk1},''')'])
                eval([MultFld1{kk1},'.',MultFld2{kk1},'=',MultFld2{kk1},'(',num2str(kk),',:);']);
            end
        end
        
        if ~isfield(Exp, 'Range') & ~isempty(hand.src.y), Exp.Range = ([hand.src.ax.x(1) hand.src.ax.x(end)] + safeget(hand.src.ax, 'dx', 0))/Opt.unit_cnv ; end
        if handles.is2D
            [B, B1, Sp] = feval(sim, Sys, Exp, Opt);
            x1 = B1'*Opt.unit_cnv;
        else
            [B, Sp] = feval(sim, Sys, Exp, Opt);
        end
        
        if strcmp(sim, 'saffron') 
            ax.x = B.';
            ax.xlabel = 'time, us';
            y = Sp.' - mean(Sp);
            pars.awin = 'gau';
            pars.awidth = 1;
            pars.ashift = 0;
            pars.aalpha = 2;
            pars.zerofill = 1;
            pars.fft = 1;
            part.opt = 'real';
            [rax, ry] = fftprocessing(ax, y, pars);
            Sp = abs(ry);
            B = rax.x;
        end
        if isfield(Exp,'blspot')
            if isfield(Sys, 'Nucs'), [Sys.I,Sys.gn] = nucdata(Sys.Nucs); end
            freq = Exp.Field/planck/1E9*Sys.gn*nmagn;
            switch Exp.blspot
            case 'none',
            case 'mims',
                Sp = Sp.*(sin((B-freq)*Exp.tau)).^2;    
            case 'davies',
            end
        end
        x = B'*Opt.unit_cnv;
        if handles.is2D, Spec(:,:,kk) = Sp;
        else, 
            Spec(:, kk) = Sp(:); % alsi 07.03.2005
        end
    end
    handles.sim{num}.x = x;
    handles.sim{num}.freq1 = safeget(Exp, 'mwFreq', 9.35)*1E9;
    handles.sim{num}.MultSize = MultSize;
    if handles.is2D, 
        handles.sim{num}.x1 = x1;
        amp = abs(max(max(Spec))-min(min(Spec)));
        amp(amp<1E-25) = 1;
        handles.sim{num}.amp = amp*ones(size(Spec,2));
    else
        if isfield(Shift, 'y')
            handles.sim{num}.x1 = Shift.y;
        else
            handles.sim{num}.x1 = [1:MultSize]';
        end
        handles.sim{num}.amp = abs(max(Spec, [], 1)-min(Spec, [], 1));
    end
    %%%%%%%%%%%% alsi 14.05.2007 %%%%%%%%%%%%%%%%%%%
    %% scale data only if neccessary %%%%%%%%%%%%%%%
        if isfield(Opt, 'isScaled')
            if ~Opt.isScaled
                handles.sim{num}.amp = 1;
            end
%         else
%             handles.sim{num}.amp = 1;
        end
    %%%%%%%%%%%% alsi 14.05.2007 %%%%%%%%%%%%%%%%%%%
    handles.sim{num}.MultSize = MultSize;
    if max(handles.sim{num}.amp)==0
        disp(['Experiment range: ', num2str(min(x)), '->', num2str(max(x))]);
        disp('All spectra have zero amplitude (Exp.Range???)');
        handles.sim{num}.y = Spec;
    else
        % protection
        handles.sim{num}.amp(handles.sim{num}.amp<1E-25) = 1;
        handles.sim{num}.y = Spec./handles.sim{num}.amp(ones(size(Spec, 1), 1), ones(MultSize, 1)); 
    end
    guidata(handles.figure1, handles);
    if nargout > 0, 
        varargout{1} = handles;
    end
    set(handles.figure1,'name','SIMPLUGIN')
catch    
    if ~isempty(lastwarn), disp(['Warning: ', lastwarn]); end
    if ~isempty(lasterr), disp(['Error: ', lasterr]); end
    set(handles.figure1,'name','SIMPLUGIN: Error')
end
% OptimiseSimulation(handles);
return

% --------------------------------------------------------------------
function plothis(handles)
handles = guidata(handles.figure1);
colors = handles.Color;

hand = guidata(handles.hh);

hand.out = {};

ax = hand.src.ax;
y = hand.src.y;

if handles.process, [y, ax] = processing(y, ax);
else
    ax.diff = '';
    ax.filt = '';
end

c = {'off', 'on'};
bsum = strcmp(get(handles.mShsumm, 'Checked'), 'on');
bdiff = strcmp(get(handles.mShdiff, 'Checked'), 'on');
bbasln = strcmp(get(handles.mSimSubBl, 'Checked'), 'on');
autobasln = strcmp(get(handles.mSimAutoBaseline, 'Checked'), 'on');
plotsim  = strcmp(get(handles.mShsim, 'Checked'), 'on');
x_sum = [];  y_sum = [];
errr = 0;

Str = get(handles.lbOpt, 'String');
OptRng = [];
for k=1:size(Str,1),
    try,eval([Str{k}(1:end) ';']);end
end

wf = safeget(Opm, 'wf', 0);
if ~isempty(y)
    exp_dat = AutoBaselineCorrection(y, strcmp(get(handles.mSimAutoBaseline, 'Checked'), 'on'));
    hand.out{1}.ax = ax;
    if handles.Larmfreq & isfield(handles, 'Larmor')
        [hand.out{1}.ax.x, hand.out{1}.y] = ShiftFreq(ax.x, exp_dat, handles.Larmor);
    else
        hand.out{1}.y = exp_dat;
    end
    sh = wf*[0:size(y,2)-1];
    hand.out{1}.y = hand.out{1}.y+sh(ones(size(y, 1), 1), :);
    hand.out{1}.ax.Color = [0, 0, 1]; % always blue
    hand.out{1}.ax.s = 0;
    hand.out{1}.xlabel = safeget(ax,'xlabel','?,');
    hand.out{1}.title = 'Src';
    hand.out{1}.ax.s = 0; % AlSi
    firstsimidx = 2;
else
    firstsimidx = 1;
end

bline = zeros(size(ax,1), 1);
sz2 = 1;
for kk = 1:handles.num
    if isempty(handles.sim{kk}.x) continue; end
    
    Shift = struct('x',0,'s',0,'Scale',1);
    if strcmp(handles.sim{kk}.stype, 'sim') | ~bbasln
        % load script 
        str = handles.sim{kk}.Script;
        for k=1:size(str,1),
            if ~isempty(strfind(str{k},'Shift.')) 
                try,eval([str{k}(3:end) ';']);end
            end
        end
        if strcmp(handles.sim{kk}.stype, 'sim')
            data_mult = exp(Opm.Amp/10);
        else
            data_mult = 1;
        end
        xTemp = handles.sim{kk}.x;
        yTemp = handles.sim{kk}.y;
        [sz1,sz2]=size(yTemp);
        if handles.Larmfreq & isfield(handles, 'Larmor')
            [xTemp, yTemp] = ShiftFreq(xTemp, yTemp, handles.Larmor); 
        end
        
        sh=[0:sz2-1]*wf;
        dx = safeget(Shift, 'dx', 0.0);
        
        Scale = safeget(Shift, 'Scale', 1).*data_mult';
        if handles.is2D
            Scale(end+1:length(handles.sim{kk}.x1),1)=Scale(1,1);
            hand.out{end}.ax.s = 0 ;                 % Shift in relative units
        else
            Scale(end+1:handles.sim{kk}.MultSize,1)=Scale(1,1);
        end
        if sz2~=length(Scale)
            adDim = [1:sz2,ones(1,length(Scale)-sz2)];
            yTemp = yTemp(:,adDim).*Scale(:,ones(sz1, 1))';
            sh    = sh(adDim);
            handles.sim{kk}.x1 = handles.sim{kk}.x1(adDim);
            [sz1,sz2]=size(yTemp);
        else
            yTemp = yTemp.*Scale(:,ones(sz1, 1))';
        end
        
        if plotsim,
            hand.out{end + 1}.ax = ax;
            fldnam = fieldnames(Shift);
            for ll=1:length(fldnam), 
                hand.out{end}.ax=setfield(hand.out{end}.ax, fldnam{ll},getfield(Shift,fldnam{ll})); 
            end
            hand.out{end}.ax.dx = dx;
            
            hand.out{end}.ax.y = handles.sim{kk}.x1;
            hand.out{end}.y = yTemp + sh(ones(sz1, 1), :);
            hand.out{end}.ax.x = xTemp;
            hand.out{end}.ax.Color = safeget(Shift, 'Color', colors(kk, :));
            hand.out{end}.ax.s = Shift.s*max(data_mult) ; % Shift in relative units
            hand.out{end}.title = ['Out', num2str(kk)];
            hand.out{end}.ax.freq1 = handles.sim{kk}.freq1;
            %% Save scripts
            for ci = 1:length(handles.sim{kk}.Script)
                hand.out{end}.Script{ci} = handles.sim{kk}.Script{ci}(3:end); 
            end
            %% Add scripts as DSC
            if strcmp(get(handles.mOptAddScripts, 'Checked'),'on')
                if ~isfield(hand.out{end},'dsc') hand.out{end}.dsc = []; end;
                for ci = 1:length(handles.sim{kk}.Script)
                    hand.out{end}.dsc = setfield(hand.out{end}.dsc, ['script',num2str(ci)],...
                        handles.sim{kk}.Script{ci}(3:end)); 
                end
            end
        end
        if (bsum | bdiff), 
            if isempty(y_sum) | isempty(x_sum), 
                y_sum = yTemp; x_sum = xTemp + dx;
            else
                for kk=1:sz2
                    spl = spline(xTemp + dx, yTemp(:,kk), x_sum);
                    y_sum(:,kk) = y_sum(:,kk) + spl(:);
                end
            end
        end
    else
        bline = bline + handles.sim{kk}.y;
    end
end

if (~errr)&(bbasln)
    hand.out{1}.y = hand.out{1}.y - bline;
end
sh=[0:sz2-1]*wf;
% Add sum
    sz = size(y_sum, 2);
if (~errr)&(bsum) & max(y_sum)~=min(y_sum)
    hand.out{end + 1}.ax = hand.out{1}.ax;
    hand.out{end}.y = y_sum+sh(ones(size(y_sum, 1), 1), :);
    hand.out{end}.ax.x = x_sum;
    hand.out{end}.ax.dx = 0;
    hand.out{end}.ax.Color = handles.SumCol; 
    hand.out{end}.ax.s = 0; 
    hand.out{end}.title = 'Sum';
    %     uses the last shift
    hand.out{end}.ax.s = Opm.dys;
end
% Add diff
if (~errr)&(bdiff) & firstsimidx==2 & max(y_sum)~=min(y_sum)
    hand.out{end + 1}.ax = hand.out{1}.ax;
    for kk=1:sz
        yy(:,kk)= exp_dat(:,kk)-spline(x_sum,y_sum(:,kk),hand.out{end}.ax.x);
    end
    hand.out{end}.y = yy + sh(ones(size(exp_dat, 1), 1), :);
    hand.out{end}.ax.Color = handles.DiffCol; 
    hand.out{end}.ax.s = 0;   
    hand.out{end}.title = 'Diff';
end
script = lscript(handles);        
guidata(handles.hh, hand);
[fpath,name,ext] = fileparts(get(handles.hh, 'FileName'));
eval([name '(''SetDataSource'', 2, hand)']);

% --------------------------------------------------------------------
function mLoad_Callback(hObject, eventdata, handles)
handles = guidata(handles.figure1);

if ischar(get(handles.popSel, 'String')), 
    str = {get(handles.popSel, 'String')};
else
    str = get(handles.popSel, 'String');
end
val = get(handles.popSel, 'Value');

switch hObject
case {handles.mLoad, handles.ptChange}
    [fname,pname] = uigetfile([handles.currdir, filesep, '*.m'],'Open Script');
    if fname == 0, return; end
    handles.currdir = pname;
    handles.sim{val}.ScriptName = [pname , filesep, fname];
    [name, script] = LoadScript([pname , filesep,  fname], handles.sim{val}.stype);
    set(handles.list, 'String', script, 'Value', 1);
    
    str{val} = [num2str(val), ': ', name];
    set(handles.popSel, 'String', str);
    handles.sim{val}.Script = script;
    handles.sim{val}.MultSize = 1;
    guidata(handles.figure1, handles);
    list_Callback(hObject, [], handles);
case handles.mLoadAll
    if exist('uigetdir')==5
        directory_name = uigetdir(handles.currdir,'Select directory');
    else
        [a, directory_name] = uiputfile([handles.currdir, filesep, 'select.dir'],'Select directory');
    end
    if ~ischar(directory_name)||isempty(directory_name), return; end
    handles.currdir = [directory_name, filesep];
    files = dir([handles.currdir, filesep, '*.m']);
    scripts = {files.name};
    
    hand = guidata(handles.hh);
    stype = 'sim';
    for kk=1:size(scripts, 2)
        [name, script] = LoadScript(fullfile(handles.currdir, scripts{kk}), 'sim');
        str{handles.num + 1} = [num2str(handles.num + 1), ': ', name];
        handles.sim{handles.num + 1}.Script = script;
        handles.sim{handles.num + 1}.x = hand.src.ax.x;
        handles.sim{handles.num + 1}.y = zeros(size(hand.src.ax.x));
        handles.sim{handles.num + 1}.x1 = 0;
        handles.sim{handles.num + 1}.MultSize = 1;
        handles.sim{handles.num + 1}.stype = stype;
        handles.Color(handles.num + 1, :) = handles.Color(mod(handles.num,7)+1, :);
        handles.num = handles.num + 1;
    end
    
    str = ChangeSimNumber(str); % changing numbers in popSel list
    set(handles.list, 'String', script, 'Value', 1);
    set(handles.popSel, 'String', str, 'Value', handles.num);
    guidata(handles.figure1, handles);
    popSel_Callback(hObject, eventdata, handles);
end
guidata(handles.figure1, handles);

% --------------------------------------------------------------------
function varargout = LoadScript(fullname, stype)
fid = fopen(fullname, 'r');
[fpath,fname,ext] = fileparts(fullname); 
if fid == -1, disp('File not found'); return; end
str = {};
while feof(fid) < 1,
    st = fgetl(fid);
    if ~isempty(st) & st(1) ~= '%',
        str{end+1} = ['l-',st];
    end
end
fclose(fid);

if nargout > 0 , 
    varargout = {fname , str.'};
end

% --------------------------------------------------------------------
function SaveScript(fullname, Scr)
fid = fopen(fullname, 'w');
if fid > 0
    for kk = 1:length(Scr)
        fprintf(fid, '%s\r\n', Scr{kk}(3:end));
    end
    fclose(fid);
else
    error(['Can''t write to file ''',fullname,'''']);
end
% --------------------------------------------------------------------
function handles = ZZ_Callback(hObject, eventdata, handles)
switch hObject
case handles.mSimAutoBaseline
    c = {'off', 'on'};
    stat = strcmp(get(handles.mSimAutoBaseline, 'Checked'), 'on');
    set(handles.mSimAutoBaseline, 'Checked', [c{~stat + 1}]);
case handles.pmOptMethod
    if get(handles.pmOptMethod, 'Value') > 5 & ...
            get(handles.butCalc, 'Value')   
        check_Callback(hObject, eventdata, handles);
        return
    end
end
handles = OptimiseAmplitude(handles);
plothis(handles);

% --------------------------------------------------------------------
function mSave_Callback(hObject, eventdata, handles)
handles = guidata(handles.figure1);
hand = guidata(handles.hh);

[pathstr,name,ext] = fileparts(get(handles.figure1, 'FileName'));
handles.currdir = safeget(handles, 'currdir', pathstr);

switch hObject
case {handles.mSave, handles.ptSave}
    num = get(handles.popSel, 'Value');
    sims = get(handles.popSel, 'String');
    if ischar(sims), string = sims;     
    else
        string = sims{num};
    end
    [fname,pname] = uiputfile([handles.currdir, filesep,  string(3:end), '.m'],'Save Script');
    if fname == 0, return; end
    handles.currdir = pname;
    [fpath,fname,ext] = fileparts([pname, fname]); 
    
    if strcmp(ext ,'') , ext = '.m'; end
    SaveScript([pname, fname, ext], get(handles.list, 'String'));
case {handles.mFileCopyClipboard, handles.ptClipboard}
    num = get(handles.popSel, 'Value');
    sims = get(handles.popSel, 'String');
    if ischar(sims), string = sims;     
    else, string = sims{num};
    end
    Scr = get(handles.list, 'String');
    out = '';
    for kk = 1:length(Scr)
        out = sprintf( '%s%s\n', out, Scr{kk}(3:end));
    end
    clipboard('copy',out);
    disp(['Simplugin: script ''',sims{num},''' copied to the clipboard.']);
case handles.mSaveAll
    if exist('uigetdir')==5
        directory_name = uigetdir(handles.currdir,'Select directory');
    else
        [a, directory_name] = uiputfile([handles.currdir, filesep, 'select.dir'],'Select directory');
    end
    if ~directory_name, return; end
    handles.currdir = directory_name;
    sims = get(handles.popSel, 'String');
    if ischar(sims), 
        cnum = 1;
    else
        cnum = 1:length(sims);
    end
    
    for num = cnum
        if ischar(sims), 
            [tmp, string] = strtokstr(sims, ': ');     
        else
            [tmp, string] = strtokstr(sims{num}, ': ');
        end
        [pa,name,ext] = fileparts(string);
        SaveScript([directory_name, filesep, name, '.m'], handles.sim{num}.Script);
        disp(['Simplugin: script ''',name,''' saved to ', directory_name]);
    end
end
guidata(handles.figure1, handles);
outstr = inimanaging([pathstr, filesep, 'kazan.ini']);
outstr.Common_For_Plugins.currDir = handles.currdir;
inimanaging([pathstr, filesep,  'kazan.ini'], outstr);

% --------------------------------------------------------------------
% --- Executes on button press in check.
function check_Callback(hObject, eventdata, handles)
[fpath,n,e] = fileparts(which('kazan'));
load([fpath, filesep,  'ico.mat']);
handles = guidata(handles.figure1);
c = {'off', 'on'};
col = {[1 0.5 0.5], [0.5, 1, 0.5]};
fieldnames = {'calcerr', 'calculate'};
stat = strcmp(get(handles.mCalc, 'Checked'), 'on');
set(handles.mCalc, 'Checked', [c{~stat + 1}]);
% set(handles.butCalc, 'Value', ~stat,'BackgroundColor', ...
%     [col{~stat + 1}]);
set(handles.butCalc, 'CData', getfield(ico, fieldnames{~stat + 1}))


if ~stat, 
    %   list_Callback(hObject, eventdata, handles)
    ePars_Callback(hObject, eventdata, handles)
end

% --------------------------------------------------------------------
% --- Executes when user attempts to close figure1.
function figure1_CloseRequestFcn(hObject, eventdata, handles)
try
    [pathstr, a, b] = fileparts(get(handles.figure1, 'FileName'));
    outstr = inimanaging([pathstr, filesep,  'kazan.ini']);
    outstr.Common_For_Plugins.currDir = handles.currdir;
    outstr.Simplugin.showSum = get(handles.mShsumm, 'Checked');
    outstr.Simplugin.showDiff = get(handles.mShdiff, 'Checked');
    outstr.Simplugin.autoBaseline = get(handles.mSimAutoBaseline, 'Checked');
    outstr.Simplugin.applSourceProc = get(handles.mSimAgr, 'Checked');
    outstr.Simplugin.loadExpFreq = get(handles.mSimApplyFreq, 'Checked');
    inimanaging([pathstr, filesep, 'kazan.ini'], outstr);
end
delete(handles.figure1);
% --------------------------------------------------------------------
function popSel_Callback(hObject, eventdata, handles)
handles = guidata(handles.figure1);
num = get(handles.popSel, 'Value');
set(handles.list, 'String', handles.sim{num}.Script, 'Value', 2);
list_Callback(hObject, eventdata, handles);
set(handles.butCol, 'BackgroundColor', handles.Color(num, :));

% --------------------------------------------------------------------
function butDel_Callback(hObject, eventdata, handles)
handles = guidata(handles.figure1);
str = get(handles.popSel, 'String');
if hObject==handles.mDeleteAll
    handles.sim = {handles.sim{1}};
    outstr = ChangeSimNumber({str{1, :}}');
    handles.num = 1;
else
    num = handles.num;
    val = get(handles.popSel, 'Value');
    if num == 1, return; end
    
    % ss = 1:num; ss = ss(ss~=val); % AlSi 07.12.04
    switch val
    case 1,
        ss = 2:num;
    case num
        ss = 1:num-1;
    otherwise
        ss = [1:val-1, val+1:num];
    end
    
    handles.sim = {handles.sim{ss}};
    handles.num = handles.num - 1;
    outstr = ChangeSimNumber({str{ss, :}}');
end
set(handles.popSel, 'String', outstr, 'Value', 1);
guidata(handles.figure1, handles);
popSel_Callback(hObject, eventdata, handles);

% --------------------------------------------------------------------
function AddSim_Callback(hObject, eventdata, handles)
handles = guidata(handles.figure1);

switch hObject
case handles.mFileBaseline
    stype = 'baseline';
    title = 'Open Baseline Script';
otherwise
    stype = 'sim';
    title = 'Open Simulation Script';
end

full = get(handles.figure1, 'FileName');
[fpath,name,ext] = fileparts(full); 
handles.currdir = safeget(handles, 'currdir', fpath);
[fname, pname] = uigetfile([handles.currdir, filesep,  '*.m'],title);
handles.currdir = pname;
if fname == 0, return; end



[name, script] = LoadScript([pname fname], stype);
set(handles.list, 'String', script, 'Value', 1);

if ischar(get(handles.popSel, 'String')),
    str = {get(handles.popSel, 'String')};
else
    str = get(handles.popSel, 'String');
end

str = ChangeSimNumber(str); % changing numbers in popSel list

str{handles.num + 1} = [num2str(handles.num + 1), ': ', name];
set(handles.popSel, 'String', str, 'Value', handles.num + 1);

hand = guidata(handles.hh);
handles.sim{handles.num + 1}.Script = script;
handles.sim{handles.num + 1}.x = hand.src.ax.x;
handles.sim{handles.num + 1}.y = zeros(size(hand.src.ax.x));
handles.sim{handles.num + 1}.x1 = 0;
handles.sim{handles.num + 1}.stype = stype;
handles.Color(handles.num + 1, :) = handles.Color(mod(handles.num,7)+1, :);
handles.num = handles.num + 1;

guidata(handles.figure1, handles);
popSel_Callback(hObject, eventdata, handles);

% --------------------------------------------------------------------
function varargout = mFileSalt_Pepper_Callback(h, eventdata, handles, varargin)
handles = guidata(handles.figure1);
func = {'Pepper', 'Salt'}; % For future, when Num. of Functions will be > 2
Name = get(h, 'Label');
for kk = 1:size(func, 2)
    if strcmp(func{kk}, Name), 
        set(h, 'Checked', 'on');
        handles.CalcType = lower(Name);
    else
        eval(['set(handles.',['mFile', func{kk}], ', ''Checked'', ''off'');']);
    end
end
guidata(h, handles);

% --------------------------------------------------------------------
function varargout = mFileCopy_Callback(h, eventdata, handles, varargin)
num = get(handles.popSel, 'Value');

script = get(handles.list, 'String');
set(handles.list, 'Value', 2);

if ischar(get(handles.popSel, 'String')),
    str = {get(handles.popSel, 'String')};
else
    str = get(handles.popSel, 'String');
end

str{handles.num + 1} = [str{num}, ' copy'];

str = ChangeSimNumber(str);

set(handles.popSel, 'String', str, 'Value', handles.num + 1);

hand = guidata(handles.hh);
handles.sim{handles.num + 1}.Script = script;
handles.sim{handles.num + 1}.x = hand.src.ax.x;
handles.sim{handles.num + 1}.y = zeros(size(hand.src.ax.x));
handles.sim{handles.num + 1}.stype = safeget(handles.sim{num}, 'stype', 'sim');
if handles.num + 1 > 7
    handles.Color(handles.num + 1, :) = [0, 0, 0];
end
handles.num = handles.num + 1;

guidata(handles.figure1, handles);
popSel_Callback(h, eventdata, handles);

% --------------------------------------------------------------------
% Change colors of components
function mSimCol_Callback(h, eventdata, handles)
handles = guidata(handles.figure1);
switch h
case handles.mColSum
    c = uisetcolor(handles.SumCol);
    if size(c, 2) > 1, handles.SumCol = c; 
    else, return;
    end
case handles.mColDiff
    c = uisetcolor(handles.DiffCol);
    if size(c, 2) > 1, handles.DiffCol = c;
    else, return;
    end        
case handles.butCol
    num = get(handles.popSel, 'Value');
    c = uisetcolor(handles.Color(num, :));
    if size(c, 2) > 1
        handles.Color(num, :) = c;
        set(handles.butCol, 'BackgroundColor', c, 'ForegroundColor', [1 1 1] - c);
    end
end
guidata(handles.figure1, handles);
plothis(handles);

% --------------------------------------------------------------------
function butCalc_Callback(hObject, eventdata, handles)
opt_type = get(handles.pmOptMethod, 'Value');
if opt_type < 7, check_Callback(hObject, eventdata, handles);
else
    set(handles.mCalc, 'Checked', 'off');
    set(handles.butCalc, 'Value', 0);
    OptimiseSimulation(handles); 
end

% --------------------------------------------------------------------
function script = lscript(handles)
handles = guidata(handles.figure1);

strings{1} = 'figure; hold on;';
for kk = 1:size(handles.sim, 2)
    for j = 1:size(handles.sim{kk}.Script, 1)
        strings{end + 1} = handles.sim{kk}.Script{j};
    end
    
    strings{end + 1} = 'if ~isfield(Exp, ''to_mT''), Exp.to_mT = 10; end';
    strings{end + 1} = 'if ~isfield(Exp, ''Range''), Exp.Range = [hand.src.ax.x(1) hand.src.ax.x(end)]/Exp.to_mT; end';
    strings{end + 1} = 'if ~isfield(Exp, ''Points'') Exp.nPoints = length(hand.src.ax.x); end';
    strings{end + 1} = 'if ~isfield(Opt, ''Verbosity''), Opt = struct(''Verbosity'', 0); end';
    
    switch handles.CalcType
    case 'Pepper'
        strings{end + 1} = '[B, Sp] = pepper(Sys, Exp, Opt);';
        strings{end + 1} = 'x = B*Exp.to_mT;';
    case 'Salt'
        strings{end + 1} = '[B, Sp] = salt(Sys, Exp, Opt);';
        strings{end + 1} = 'x = B;';
    end
    %   strings{end + 1} = ['Scale = ', num2str(handles.Scale(kk)), ';'];
    strings{end + 1} = ['Zoom = ', num2str(exp(handles.div/10)), ';'];
    strings{end + 1} = 'y = Sp/abs(max(Sp)-min(Sp))*Scale*Zoom;';
    strings{end + 1} = ['Color = ', num2str(handles.Color(kk, :)), ';'];
    %   strings{end + 1} = ['Shift = ', num2str(handles.Shift(kk)), '*Zoom;']; % Shift in relative units
    strings{end + 1} = 'plot(x, y+Shift.y, ''Color'', Color);';
    
end

strings{end + 1} = 'hold off;';

script = strings';

% --------------------------------------------------------------------
function outstr = ChangeSimNumber(str)
for k = 1:size(str, 1)
    tempstr = str{k};
    if length(tempstr)>2,
        [strstart, strend] = strtokstr(tempstr, ': ');
    else
        strend = tempstr;
    end
    if ~isempty(strend), str{k} = [num2str(k), ': ', strend];
    else
        str{k} = [num2str(k), ': ', strstart];
    end
end
outstr = str;

% --------------------------------------------------------------------
function mSimChi_Callback(h, eventdata, handles)
hand = guidata(handles.hh);
ax = hand.src.ax;
y = hand.src.y;

if handles.process, [y, ax] = processing(y, ax);end

Str = get(handles.lbOpt, 'String');
for k=1:size(Str,1),
    try,eval([Str{k}(3:end) ';']);end
end

y_sum = 0;
for k = 1:length(handles.sim) 
    str = handles.sim{k}.Script;
    for m=1:size(str,1),
        if ~isempty(strfind(str{m},'Shift.')) 
            try,eval([str{m}(3:end) ';']);
            catch, disp(lasterror);
            end
        else, Shift.Scale = 1;
        end
    end
    if strcmp(handles.sim{k}.stype, 'sim')
        data_mult = exp(Opt.Amp/10);
    else
        data_mult = 1;
    end
    y_sum = handles.sim{k}.y*Shift.Scale*data_mult + y_sum;
end

diff = y - y_sum;

set(hand.stCoordinates,'String', ['Chi^2: ', num2str(sum(diff.^2))]);

% --------------------------------------------------------------------
function mSimAgr_Callback(h, eventdata, handles)
str = {'off', 'on'};
stat = strcmp(get(handles.mSimAgr, 'Checked'), 'on');
set(handles.mSimAgr, 'Checked', str{~stat + 1});
handles.process = ~stat;
guidata(handles.figure1, handles);
plothis(handles);

% --------------------------------------------------------------------
function mFileAddPar_Callback(h, eventdata, handles)
handles = guidata(handles.figure1);
scr = get(handles.list, 'String');
pos = get(handles.list, 'Value');

num = get(handles.popSel, 'Value');
answer = kv_constructfield;

if isempty(answer) , return; end

scrlen = length(scr);
script = {};
for count = 1:pos
    script{end+1} = scr{count};
end
script{end+1} = ['l-',answer];
if pos ~= scrlen,
    for count = pos+1:scrlen
        script{end+1} = scr{count};
    end
end

set(handles.list, 'String', script, 'Value', pos+1);
handles.sim{num}.Script = script;
guidata(handles.figure1, handles);
list_Callback(h, [], handles);

% --------------------------------------------------------------------
function mFileChPar_Callback(h, eventdata, handles)
handles = guidata(handles.figure1);
scr = get(handles.list, 'String');
pos = get(handles.list, 'Value');

num = get(handles.popSel, 'Value');

answer = kv_constructfield(scr{pos}(3:end));
if isempty(answer) , return; end

scr{pos} = [scr{pos}(1:2),answer];

set(handles.list, 'String', scr, 'Value', pos);
handles.sim{num}.Script = scr;
guidata(handles.figure1, handles);
list_Callback(h, [], handles);

% --------------------------------------------------------------------
function mFileDelPar_Callback(h, eventdata, handles)
handles = guidata(handles.figure1);
scr = get(handles.list, 'String');
pos = get(handles.list, 'Value');

num = get(handles.popSel, 'Value');

button = questdlg('Do you want to continue?',...
    'Delete parameter','Yes','No','Yes');
if strcmp(button,'No'), return; end

scrlen = length(scr);
script = {};
for count = 1:pos-1,
    script{end + 1} = scr{count};
end
if pos~=scrlen,
    for count = pos+1:scrlen
        script{end + 1} = scr{count};        
    end
end

set(handles.list, 'String', script, 'Value', pos);
handles.sim{num}.Script = script;
guidata(handles.figure1, handles);

% --------------------------------------------------------------------
function mFileChName_Callback(h, eventdata, handles)
handles = guidata(handles.figure1);
scr = get(handles.popSel, 'String');
if ischar(scr), scr = {scr}; end

num = get(handles.popSel, 'Value');

prompt = {'Type name'};
dlg_title = 'Change simulation name';
num_lines= 1;
if length(scr{num}) > 2, 
    [en, st] = strtokstr(scr{num}, ': ');
else
    st = scr{num};
end

def     = {st};
answer  = inputdlg(prompt,dlg_title,num_lines,def);

scr{num} = [num2str(num), ': ',answer{1}];
set(handles.popSel, 'String', scr, 'Value', num);

% --------------------------------------------------------------------
function mAbout_Callback(hObject, eventdata, handles)
handles = guidata(handles.figure1);

if hObject==handles.mAboutHelp
    [fpath,n,e] = fileparts(which('kazan'));
    filename = ['file:///',fpath, filesep, 'help', filesep, 'simplugin_sims.html'];
    filename(filename=='\')='/';
    stat = web(filename);
else
    msgbox({'simplugin';
        'EPR/ENDOR simulations interface';
        'by Boris Epel and Alexey Silakov, 2003-04';
        'EasySpin (www.esr.ethz.ch) syntax rules are used';'';
        'Syntax:';
        'Sys/Opt/Exp/Shift.fld_name = any_Matlab_expression';
        'i.e. Opt.ThetaRange([1,2])=[0,pi]';'';
        'Parameters additional to EasySpin:';
        'Opt.Sim - any sim of the type [x,y]=sim_name(Sys, Exp, Opt)';
        '  = pepper, salt, sugar, sugar1 <pepper>';
        '''2D'' in the name of simulation result in [x,x1,y]=sim_name(Sys, Exp, Opt)';
        '  = triplesugar2D';
        'Exp.gField - Field in g-factor units (salt, sugar)';'';
        'Shift parameters:';
        's - simulation shift along x-axis <0>';
        'Scale - amplitude of simulation <1>';
        'wf - waterfall (shift of 2D slices) <1>';
        'line <''-''> and color [r g b]';'';
        'Some parameters (see line 314) could be arrays (COLUMNS!!!)';
        'i.e.: Sys.gField=[2.001;2.002;2.003] or Opt.ThetaRange(:,1) = [0:8]''*pi/15';
        'Shift.y - axes';
        'Shift.yslicelabel - label format (''%5.2f'')';
        'Shift.yslicelabeldx';
        'Shift.yslicelabeldy '; '';
        'Opt.unit_cnv - convert units, e.g. G to mT and back. for "pepper" it is done automatically';
        'Opt.isScaled = 0, disables amplitude normalisation';
        
    }, ...
        'About', 'help')
end
% --------------------------------------------------------------------
function varargout = mSimApplyFreq_Callback(h, eventdata, handles, varargin)
str = {'off', 'on'};
stat = strcmp(get(handles.mSimApplyFreq, 'Checked'), 'on');
set(handles.mSimApplyFreq, 'Checked', str{~stat + 1});
handles.applyfreq = ~stat;
guidata(handles.figure1, handles);

% --------------------------------------------------------------------
function Str1 = UpdatePars(Str, par, val)
if ischar(Str), Str={Str}; end;
Str2=[];
if ischar(val)
    Str2{end+1}=[par,' = ',val];
else
    if length(val)==1, Str2{end+1}=[par,' = ',num2str(val)]; 
    else
        for kk=1:size(val,1)
            for ll=1:size(val,2)
                if size(val,1)==1
                    Str2{end+1}=[par,'(',num2str(ll),') = ',num2str(val(kk,ll))]; 
                else
                    Str2{end+1}=[par,'(',num2str(kk),', ',num2str(ll),') = ',num2str(val(kk,ll))]; 
                end
            end
        end
    end
end
max_insert = length(Str2);
cur_ins=1;

Str1=[];
done  = 0;
for kk=1:length(Str)
    if isempty(strfind(Str{kk},par))
        Str1{end+1}=Str{kk};
    else
        if done, continue; end
        done=1;
        stat = Str{kk}(1:2);
        if cur_ins<=max_insert
            if ~(strcmp(stat(1), 'l')|strcmp(stat(1), 'c'))| ~(strcmp(stat(2), '-')|strcmp(stat(1), 'o'))
                Str1{end+1} = Str2{cur_ins};
            else
                Str1{end+1}=[stat,Str2{cur_ins}]; %Str1{end+1}=[Str2{cur_ins}]; %
            end
            cur_ins=cur_ins+1;
        end
        
    end
end

for kk=cur_ins:max_insert, Str1{end+1}=[Str2{kk}]; end % 'l-'

% --------------------------------------------------------------------
function varargout = mRefresh_Callback(h, eventdata, handles, varargin)
plothis(handles);

% --------------------------------------------------------------------
function varargout = mSortPars_Callback(h, eventdata, handles, varargin)
scr = sort(get(handles.list, 'String'));
num = get(handles.popSel, 'Value');

set(handles.list, 'String', scr, 'Value', 1);
handles.sim{num}.Script = scr;
guidata(handles.figure1, handles);
list_Callback(h, [], handles);


% --------------------------------------------------------------------
function varargout = pbOptimise_Callback(h, eventdata, handles, varargin)

% --------------------------------------------------------------------
function varargout = popup_Callback(h, eventdata, handles, varargin)
num = get(handles.popup, 'Value');
str = get(handles.popup, 'String');
dim = str2num([str{num}]);
set(handles.ePars, 'Step', dim);

% --------------------------------------------------------------------
function varargout = slOpt_Callback(h, eventdata, handles, varargin)
shift = get(handles.slOpt, 'Value');
set(handles.slOpt, 'Value', 0.5);
num = get(handles.pmOpt, 'Value');
str = get(handles.pmOpt, 'String');
dim = str2num([str{num}]);

CurVal = str2num(get(handles.eOpt, 'String'));
Val = CurVal + dim*(shift - 0.5)*2;

set(handles.eOpt, 'String', num2str(Val, 15));
eOpt_Callback(h, eventdata, handles, varargin);

% --------------------------------------------------------------------
function varargout = pmOpt_Callback(h, eventdata, handles, varargin)
num = get(handles.pmOpt, 'Value');
str = get(handles.pmOpt, 'String');
dim = str2num([str{num}]);
set(handles.eOpt, 'Step', dim);

% --------------------------------------------------------------------
function varargout = eOpt_Callback(h, eventdata, handles, varargin)
handles = guidata(handles.figure1);
Res = get(handles.lbOpt, 'String');
val = get(handles.lbOpt, 'Value');
newval = trim(get(handles.eOpt, 'String'));
[name, str1] = strtok(Res{val}, '=');
Res{val}=[name,'= ', num2str(newval)];
set(handles.lbOpt, 'String', Res);

if strfind(name, 'OptRng'), SetRange(handles); 
end
plothis(handles);

% --------------------------------------------------------------------
function varargout = lbOpt_Callback(h, eventdata, handles, varargin)
Res = get(handles.lbOpt, 'String');
val = get(handles.lbOpt, 'Value');
[name, str1] = strtok(Res{val}, '=');
str2 = [trim(str1(3:end)), ' '];
set(handles.eOpt, 'String', str2);

% --------------------------------------------------------------------
function handles = OptimiseAmplitude(handles)
opt_type = get(handles.pmOptMethod, 'Value');
if opt_type < 2 | opt_type > 6, return;
end

str = get(handles.lbOpt, 'String');
OptRng=[];
for kk=1:size(str,1),
    try,eval([str{kk}(3:end) ';']);end
end

try
simnum = size(handles.sim,2);
hand = guidata(handles.hh);
if simnum > 0 & ~isempty(hand.src.y)
    idx = handles.DataRange.idx;
    y = real(AutoBaselineCorrection(hand.src.y, strcmp(get(handles.mSimAutoBaseline, 'Checked'), 'on')));
    switch opt_type
    case 2, % Optimise spread
        for kk=1:length(handles.sim)
%             max1 = max(real(handles.sim{kk}.y(handles.sim{kk}.idx,:)));
%             amp1(kk,:) = abs(max1-min(real(handles.sim{kk}.y(handles.sim{kk}.idx,:))));
            max1 = max(real(handles.sim{kk}.y));
            amp1(kk,:) = abs(max1-min(real(handles.sim{kk}.y)));
        end
        max2 = max(y(idx,:));
        amp2 = abs(max2-min(y(idx,:)));
    case 3, % Optimise maximum
        for kk=1:length(handles.sim)
%             amp1(kk,:) = max(real(handles.sim{kk}.y(handles.sim{kk}.idx,:)));
            amp1(kk,:) = max(real(handles.sim{kk}.y));
        end
        amp2 = max(y(idx,:));
    case 4, % Optimise integral
        for kk=1:length(handles.sim)
            amp1(kk,:) = sum(real(handles.sim{kk}.y(handles.sim{kk}.idx,:)));
            amp1(kk,:) = amp1(kk,:) * diff(handles.sim{kk}.x(1:2));
        end
        amp2 = sum(y(idx,:)) * diff(hand.src.ax.x(1:2));
    case 5, % Optimize doubleintegral
        for kk=1:length(handles.sim)
            amp1(kk,:) = sum(cumsum(real(handles.sim{kk}.y(handles.sim{kk}.idx,:))));
            amp1(kk,:) = amp1(kk,:) * (diff(handles.sim{kk}.x(1:2)))^2;
        end
        amp2 = sum(cumsum(y(idx,:))) * (diff(hand.src.ax.x(1:2)))^2;
    case 6, % Preserve double integral
        for kk=1:length(handles.sim)
            amp1(kk,:) = sum(cumsum(real(handles.sim{kk}.y(handles.sim{kk}.idx,:))));
            amp1(kk,:) = amp1(kk,:) * (diff(handles.sim{kk}.x(1:2)))^2;
        end
        amp2 = sum(cumsum(y(idx,:))) * (diff(hand.src.ax.x(1:2)))^2;
    otherwise 
        amp1=1; amp2=1;
    end
    if isempty(amp1), amp1=amp2; end;
    
    if ~sum(~(amp1(1,:)))
        res = 10*log(amp2./amp1(1,:));
    else
        res  =1;
    end
    if isempty(res), res=1; end
    
    for kk=2:length(handles.sim)
        cc = amp1(1,:)./amp1(kk,:);
        handles.sim{kk}.y = handles.sim{kk}.y.*cc(ones(size(handles.sim{kk}.y,1),1),:);
    end
else
    res = 0;
end

if(opt_type~=6)
    Str = get(handles.lbOpt, 'String');
    val = min(length(Str), get(handles.lbOpt, 'Value'));
    set(handles.lbOpt, 'String', UpdatePars(Str, 'Opm.Amp', res),'Value',val);
end
catch
    error(sprintf('%d simulations and %d experiments', length(handles.sim), length(amp1)))
end
guidata(handles.figure1, handles)

% --------------------------------------------------------------------
function y = AutoBaselineCorrection(y, mode)
switch mode
case 1,
    if size(y,1) > 15
        bsl = sum(y([1:5, end-4:end], :), 1)/10;
        y = y - bsl(ones(size(y,1),1),:);
    end    
end

% --------------------------------------------------------------------
function mContOptSet_Callback(h,eventdata,handles)
Res = get(handles.list, 'String');
val = get(handles.list, 'Value');

switch h
case handles.mContOpt, Res{val}(2)='+';
case handles.mContNotOpt, Res{val}(2)='-';
end
set(handles.list, 'String', Res);
handles.sim{get(handles.popSel, 'Value')}.Script = Res;
guidata(handles.figure1, handles)


% --------------------------------------------------------------------
function varargout = mContext1_Callback(h, eventdata, handles, varargin)
% menu recreating
delete(findobj(h, 'UserData', 'ranges'));
if ~isempty(handles.DataRange.bracket)
    for k=1:size(handles.DataRange.bracket,1)
        uimenu(h, 'Label', ['Delete Range ', num2str(k)], 'UserData', 'ranges',...
            'Tag', 'mPluginsReposition', ...
            'Callback', ['simplugin(''mCont1DelRange_Callback'',gcbo,',num2str(k),',guidata(gcbo))']);
    end
end

% --------------------------------------------------------------------
function varargout = mCont1AddRange_Callback(h, eventdata, handles, varargin)
str = get(handles.lbOpt, 'String');
hand = guidata(handles.hh);
if isfield(hand.src, 'ax'),
    numends  = num2str(max(hand.src.ax.x)); 
    numbegs = num2str(min(hand.src.ax.x));
else
    numends  = '1'; numbegs = '1';
end
str{end + 1} = ['**OptRng(end+1, 1) = ',numbegs];
str{end + 1} = ['**OptRng(end, 2) = ', numends];
guidata(handles.figure1, handles);
set(handles.lbOpt, 'String', str);

SetRange(handles);
plothis(handles);

% --------------------------------------------------------------------
function varargout = mCont1DelRange_Callback(h, del, handles, varargin)
str = get(handles.lbOpt, 'String');
idx = 1:size(handles.DataRange.bracket,1);
idx = idx(idx~=del);
handles.DataRange.bracket = handles.DataRange.bracket(idx,:);
str = UpdatePars(str, 'OptRng', handles.DataRange.bracket);
set(handles.lbOpt, 'String',str, 'Value', 1);
guidata(handles.figure1, handles);

ZZ_Callback(h, [], handles);

% --------------------------------------------------------------------
function handles = SetRange(handles)
str = get(handles.lbOpt, 'String');
OptRng = [];
for kk=1:length(str)
    eval([str{kk}(3:end),';']);
end
handles.DataRange.bracket = OptRng;

handles.DataRange.idx = [];
hand = guidata(handles.hh);
if ~isempty(hand.src.ax.x)
    x=hand.src.ax.x;
    if isempty(OptRng), handles.DataRange.idx = ones(length(x),1) > 0;
    else, handles.DataRange.idx = zeros(length(x),1) > 0;
    end
    for kk=1:size(OptRng,1)
        handles.DataRange.idx = (x >= OptRng(kk,1) &  x <= OptRng(kk,2)) ...
            | handles.DataRange.idx;
    end
end
for kk=1:length(handles.sim)
    x=handles.sim{kk}.x;
    if isempty(OptRng), handles.sim{kk}.idx = ones(length(x),1) > 0;
    else, handles.sim{kk}.idx = zeros(length(x),1) > 0;
    end
    for ll=1:size(OptRng,1)
        handles.sim{kk}.idx = (x >= OptRng(ll,1) &  x <= OptRng(ll,2)) ...
            | handles.sim{kk}.idx;
    end
end
handles = OptimiseAmplitude(handles);
guidata(handles.figure1, handles);

% --------------------------------------------------------------------
function handles = OptimiseSimulation(handles)
opt_type = get(handles.pmOptMethod, 'Value');
if opt_type < 6, return;
end
try
    set(handles.figure1,'name','Busy')
    hand = guidata(handles.hh);
    nsim = length(handles.sim);
    
    str = get(handles.lbOpt, 'String');
    OptRng=[];
    for kk=1:size(str,1),
        try,eval([str{kk}(3:end) ';']);end
    end
    
    BigOpm.Amp = exp(Opm.Amp/10);
    BigOpm.Data = AutoBaselineCorrection(real(hand.src.y), strcmp(get(handles.mSimAutoBaseline, 'Checked'), 'on'));
    BigOpm.DataIdx = handles.DataRange.idx;
    
    % amplitude arguments 
    nsimargs = nsim;
    simargs(1:nsim)=1;
    
    % construction of  spin system
    for kk=1:nsim
        % load script 
        Sys = []; Exp = []; Opt = []; Shift = [];
        str = handles.sim{kk}.Script;
        BigOpm.Vars{kk}=[];
        for k=1:length(str),
            if str{k}(2)=='+'
                nsimargs = nsimargs + 1;
                [key, str1] = strtok(str{k}, '=');
                val = [trim(str1(3:end)), ' '];
                simargs(end+1)=str2num(val);
                BigOpm.Vars{kk}{end+1}=[key(3:end),'=x(',num2str(nsimargs),');'];
            else
                try,eval([str{k}(3:end) ';']);
                catch, disp([str{k}, ': wrong line']);    
                end
            end
        end
        
        simargs(kk) = simargs(kk) * Shift.Scale;
        % Default values 
        ll = length(hand.src.ax.x);
        if ll < 5, ll=1024; end; 
        Exp.nPoints = safeget(Exp, 'nPoints', ll);
        Opt.Verbosity = safeget(Opt, 'Verbosity', 0);
        
        sim = safeget(Opt, 'Sim', 'pepper');
        
        if safeget(Exp, 'gField', 0) > .5
            Exp.Field = g2fld(Exp.gField, Exp.mwFreq*1E9)*1E3;
        end 
        
        MultFld1 = {'Exp', 'Exp', 'Opt', 'Exp','Sys', 'Shift'};
        MultFld2 = {'Field','gField', 'ThetaRange', 'ExciteWidth','lwEndor', 'Scale'};
        MultSize = 1;
        % Determine multiple argument variables 
        for jj=1:size(MultFld1, 2)
            if eval(['isfield(',MultFld1{jj},',''',MultFld2{jj},''')'])
                eval([MultFld2{jj},'=',MultFld1{jj},'.',MultFld2{jj},';']);
                MultSize = eval(['max([MultSize, size(',MultFld2{jj},',1)]);']);
            end
        end
        
        % Construct correct-dimensioned variables
        for jj=1:size(MultFld1, 2)
            if eval(['exist(''',MultFld2{jj},''')'])==1
                eval(['sz = size(',MultFld2{jj},',1);']);
                if sz < MultSize
                    eval([MultFld2{jj},'(',num2str(sz+1),':',num2str(MultSize),',:)=',...
                            MultFld2{jj},'(zeros(',num2str(MultSize-sz),',1)+end,:);']);
                end
            end
        end
        
        is2D = ~isempty(strfind(sim, '2D'));
        for jj=1:MultSize
            % Determine multiple argument variables 
            for jj1=1:size(MultFld1, 2)
                if eval(['isfield(',MultFld1{jj1},',''',MultFld2{jj1},''')'])
                    eval([MultFld1{jj1},'.',MultFld2{jj1},'=',MultFld2{jj1},'(',num2str(jj),',:);']);
                end
            end
            
            if strcmp(sim, 'pepper'), defltunit = 10; else defltunit = 1; end 
            Opt.unit_cnv = safeget(Opt, 'unit_cnv', defltunit); % mt -> G conversion
            if ~isfield(Exp, 'Range') & ~isempty(hand.src.y), Exp.Range = ([hand.src.ax.x(1) hand.src.ax.x(end)] + hand.src.ax.dx)/Opt.unit_cnv ; end
            
            BigSys{jj,kk}=Sys;
            BigExp{jj,kk}=Exp;
            BigOpt{jj,kk}=Opt;
        end
    end
catch
    disp(lasterr)
end

switch opt_type
case 8, optset = optimset('MaxIter',10);
case 9, optset = optimset('MaxIter',100);
case 10, optset = optimset('TolFun',1e-6);
case 11, optset = optimset('TolFun',1e-8);
otherwise, optset = [];
end

% newsimargs = fminsearch(@Opt_function, simargs, optset, BigSys, BigExp, BigOpt, BigOpm);
err = .2;
newsimargs = fmincon(@Opt_function, simargs, ...
    [], [], [], [], simargs*(1-err), simargs*(1+err), nonlcon, ...
optset, BigSys, BigExp, BigOpt, BigOpm);
disp(['Shi:',num2str(Opt_function(newsimargs, BigSys, BigExp, BigOpt, BigOpm))])
disp(simargs)
disp(newsimargs)

set(handles.figure1,'name','simplugin');
guidata(handles.figure1, handles);

% --------------------------------------------------------------------
function [c,ceq] = nonlcon(x)
c = [];
ceq = [];

% --------------------------------------------------------------------
function shi=Opt_function(x, Sy, Ex, Op, Opm)
% x structure
% 1..nsimulations     amplitudes of sim
% nsimulations+1..    fitting parameters

% Opm.Vars{}={....}   parameters needed 
%                     for particular simulation
% Opm.Data     experimental data 
% Opm.DataIdx  index of points
try
sims = size(Sy,2);
% slice cycle 
shi = 0;
for kk=1:size(Sy,1)
    SSp=0;
    % sim cycle
    for ll=1:sims
        Sys = Sy{kk,ll};
        Exp = Ex{kk,ll};
        Opt = Op{kk,ll}; 
        % set variable parameters
        for jj=1:size(Opm.Vars{ll},2)
            eval(Opm.Vars{ll}{jj})
        end
        try, 
            [B, Sp] = feval(Opt.Sim, Sys, Exp, Opt);
            SSp=SSp+Sp*x(ll);
        catch, SSp = 1E56;
        end
    end
    shi = shi + sum((Opm.Data(Opm.DataIdx,kk)-SSp(Opm.DataIdx)'*Opm.Amp(kk)).^2);
end
catch, shi = 1E65;
end
% --------------------------------------------------------------------
function mFileMove_Callback(hObject, eventdata, handles)
handles = guidata(hObject);
Str = get(handles.list, 'String');
Val = get(handles.list, 'Value');
if ischar(Str), return; end % ischar==1 means only 1 string in the list.
maxval = length(Str);
Str1 = Str;
switch hObject
    case {handles.mFileMoveUp, handles.mMoveUp}
        if Val == 1, return; end
        Str1{Val-1} = Str{Val};
        Str1{Val} = Str{Val-1};        
        newval = Val-1;
    case {handles.mFileMoveDown, handles.mMoveDown}
        if Val == maxval, return; end
        Str1{Val+1} = Str{Val};
        Str1{Val} = Str{Val+1};                
        newval = Val+1;
end
set(handles.list, 'String', Str1, 'Value', newval);
% AAAAAAAAAAA no update 

% % --------------------------------------------------------------------
function varargout = figure1_ResizeFcn(h, eventdata, handles, varargin)
% % 140 210
FigaPos = get(handles.figure1, 'Position');

if FigaPos(3)<130, FigaPos(3) = 130; set(handles.figure1, 'Position', FigaPos); end
if FigaPos(4)<210, FigaPos(4) = 210; set(handles.figure1, 'Position', FigaPos); end
brd = 7; % board
butCalcH = 30; % points
UpperFrameH = 25;
EditH = 25;
ButLR= 20;
TopFigure = 20+brd;
CheckH = 25;
Calc2Pos = [brd, brd, butCalcH, butCalcH]; % Calculate
set(handles.butCalc,'Units', 'pixels', 'Position', Calc2Pos);

CalcPos = [2*brd + butCalcH, brd, butCalcH, butCalcH]; % Calculate All

set(handles.pbCalcAll, 'Units', 'pixels', 'Position', CalcPos);
ReCalcPos = [3*brd + 2*butCalcH, brd, butCalcH, butCalcH]; % Calculate
set(handles.pbReCalc,'Units', 'pixels', 'Position', ReCalcPos);

%  Optimization
BotFramePos = [brd, sum(Calc2Pos([2, 4]))+brd, FigaPos(3)-2*brd, 120]; 
set(handles.frame2,'Units', 'pixels', 'Position', BotFramePos);

SlUD_x= 17; SlUD_y = EditH;
Pop2Pos = [BotFramePos(1) + brd, BotFramePos(2)+brd+2, ...
        (BotFramePos(3) - 4*brd - SlUD_x)/2, EditH];
set(handles.pmOpt,'Units', 'pixels', 'Position', Pop2Pos);

eResPos = [Pop2Pos(1) + brd + Pop2Pos(3), BotFramePos(2)+brd+2, ...
         (BotFramePos(3) - 4*brd - SlUD_x)/2, EditH];
set(handles.eOpt,'Position', eResPos, 'Units', 'pixels');


% set(handles.slOpt,'Units', 'pixels', 'Position', SlPos2);

lbOptPos = [Pop2Pos(1), sum(Pop2Pos([2, 4])) + brd, BotFramePos(3)-2*brd, BotFramePos(4) - 45-4*brd];
set(handles.lbOpt,'Units', 'pixels', 'Position', lbOptPos);

pmMethPos = [Pop2Pos(1), sum(lbOptPos([2, 4])) + brd, BotFramePos(3)-2*brd, 15];
set(handles.pmOptMethod,'Units', 'pixels', 'Position', pmMethPos);

textOptPos = [Pop2Pos(1) + brd, sum(BotFramePos([2, 4]))-brd, 50, 14];
set(handles.text3,'Units', 'pixels', 'Position', textOptPos);

% Simulation

FramePos = [brd, brd + BotFramePos(2) + BotFramePos(4),...
        FigaPos(3)-2*brd, FigaPos(4)-sum(BotFramePos([2,4]))-3*brd];
set(handles.frame1,'Units', 'pixels', 'Position', FramePos);

EditPos = [FramePos(1)+brd, FramePos(2)+brd,FramePos(3)-2*brd, EditH];
set(handles.ePars,'Position', EditPos, 'Units', 'pixels');

SlPos = [EditPos(1), EditPos(2) + EditPos(4)+brd, 25, 25];
% set(handles.slPars, 'Units', 'pixels', 'Position', SlPos);

PopPos = [SlPos(1) + SlPos(3)+brd, EditPos(2)+EditPos(4)+brd-2,...
         EditPos(3)-brd - SlPos(3),EditH];
set(handles.popup, 'Units', 'pixels', 'Position', PopPos);

ListPos = [SlPos(1), SlPos(2)+SlPos(4)+brd, EditPos(3),...
         max(FramePos(4)-EditPos(4)-SlPos(4)-18-5*brd, 1)];
set(handles.list, 'Units', 'pixels', 'Position', ListPos);

PopSelPos = [FramePos(1)+brd, sum(ListPos([2, 4]))+brd-2, ...
         3*(FramePos(3)-2*brd)/4-brd, EditH];
set(handles.popSel, 'Units', 'pixels', 'Position', PopSelPos);

ColPos = [FramePos(1)+ PopSelPos(3) + 2*brd, sum(ListPos([2, 4]))+brd, ...
         (FramePos(3)-2*brd)/4, EditH]; 
set(handles.butCol, 'Units', 'pixels', 'Position', ColPos);

textSimPos = [PopSelPos(1) + brd, sum(FramePos([2, 4])), 50, 14];
set(handles.text4,'Units', 'pixels', 'Position', textSimPos);

% --------------------------------------------------------------------
function commout = findcom(handles)
% finding common strings in several scripts
handles = guidata(handles.figure1);
num = length(handles.sim);
refscr = handles.sim{1}.Script;
for cc = 1:num-1 % all simulations
    counter = [];
    for ci = 1:length(refscr) % all strings in reference list
        [r_names, r_val] = strtok(refscr{ci}(3:end), '=');
        [r_str, r_name] = strtok(r_names, '.');
        for cj = 1:length(handles.sim{cc+1}.Script);  % all strings in current list
            [c_names, c_val] = strtok(handles.sim{cc+1}.Script{cj}(3:end), '=');
            [c_str, c_name] = strtok(c_names, '.');
            if strcmp(r_str, c_str)&strcmp(trim(r_name(2:end)), trim(c_name(2:end)))&strcmp(trim(r_val(2:end)), trim(c_val(2:end)))
               counter(end+1) = ci;
            end
        end
    end
    if ~isempty(counter)
        refscr = refscr(counter);
    else
        refscr = {};
        break;
    end
end
commout = refscr;

% --------------------------------------------------------------------
function handles = change_comm(handles, commstr)
num = length(handles.sim);
[r_names, r_val] = strtok(commstr(3:end), '=');
[r_str, r_name] = strtok(r_names, '.');
if num == 1, disp('change_comm: Only one script is loaded.'); return; end
for cc = 1:num % all simulations
    counter = 0;
    for cj = 1:length(handles.sim{cc}.Script);  % all strings in current list
        [c_names, c_val] = strtok(handles.sim{cc}.Script{cj}(3:end), '=');
        [c_str, c_name] = strtok(c_names, '.');
        if strcmp(trim(r_str), trim(c_str))&strcmp(trim(r_name(2:end)), trim(c_name(2:end)))
            counter = cj;
            break;
        end
    end
    if strcmp(commstr(1), 'c'),
        if counter,
            handles.sim{cc}.Script{cj} = commstr;
        else
            handles.sim{cc}.Script{end+1} = commstr;
        end
    else
        if counter,
            handles.sim{cc}.Script{cj} = ['l', handles.sim{cc}.Script{cj}(2:end)];
        else
        end
    end
end    
% --------------------------------------------------------------------
function mConComm_Callback(h, sss, handles)
str = get(handles.list, 'String');
val = get(handles.list, 'Value');
switch h
    case handles.mConMakeCom
        handles = change_comm(handles, ['c', str{val}(2:end)]);
    case handles.mConMakeLoc
        handles = change_comm(handles, ['l', str{val}(2:end)]);
end
guidata(handles.figure1, handles);
popSel_Callback(handles.popSel, [], handles);

% --------------------------------------------------------------------
function handles = calc_all(handles)
if strcmp(get(handles.mCalc, 'Checked'), 'on'), 
for nn = 1:length(handles.sim)
    switch handles.sim{nn}.stype
        case 'sim'
            handles = pepcalc(nn, handles);
            %handles = SetRange(handles);
        case 'baseline'
            handles = bscalc(get(handles.popSel, 'Value'), handles);
    end
end
guidata(handles.figure1, handles);
plothis(handles);
end

% --------------------------------------------------------------------
function varargout = mSimLoadFreq_Callback(h, eventdata, handles, varargin)
hand = guidata(handles.hh);

% Loading of the experiment mw frequncy 
newval = safeget(hand.src.ax, 'freq1', 9.34)*1E-9;
Str = get(handles.list,'String');
Str = UpdatePars(Str, 'Exp.mwFreq', newval);
set(handles.list,'String', Str);

for kk=1:size(handles.sim,2)
    handles.sim{kk}.Script = UpdatePars(handles.sim{kk}.Script, 'Exp.mwFreq', newval);
end
guidata(handles.figure1, handles);
list_Callback(handles.figure1, [], handles);

% --------------------------------------------------------------------
function varargout = mSimLoadRange_Callback(h, eventdata, handles, varargin)
hand = guidata(handles.hh);
range = [min(hand.src.ax.x), max(hand.src.ax.x)];


% Loading of the experiment range 
Str = get(handles.list,'String');
Str = UpdatePars(Str, 'Exp.Range', range);
set(handles.list,'String', Str);

for kk=1:size(handles.sim,2)
    handles.sim{kk}.Script = UpdatePars(handles.sim{kk}.Script, 'Exp.Range', range);
end
guidata(handles.figure1, handles);
list_Callback(handles.figure1, [], handles);

% --------------------------------------------------------------------
function res = nfreq2fld(nfreq, nuc)

if nargin < 2, nuc = '1H'; end

[I, gn] = nucdata(nuc);
res = nfreq * planck / gn / nmagn;

% --------------------------------------------------------------------
function res = fld2nfreq(fld, nuc)

if nargin < 2, nuc = '1H'; end

[I, gn] = nucdata(nuc);
gn = gn(1);
res = fld*gn*nmagn/planck;


% --------------------------------------------------------------------
function varargout = mSimLarm_Callback(h, eventdata, handles, varargin)
str = {'off', 'on'};
stat = strcmp(get(handles.mSimLarm, 'Checked'), 'on');
set(handles.mSimLarm, 'Checked', str{~stat + 1});
handles.Larmfreq = ~stat;
guidata(handles.figure1, handles);
plothis(handles);

% --------------------------------------------------------------------
function [newx, newy] = ShiftFreq(x, y, shift) 

if length(shift)==size(y,2)
    maxx = max(x);
    minx = min(x);
    nminx = minx-max(shift);
    nmaxx = maxx-min(shift);
    newx = linspace(nminx,nmaxx,length(x));
    for ii=1:size(y,2)
        newy(:,ii)=spline(x-shift(ii),y(:,ii),newx)';
        idx = newx < minx-shift(ii) | newx > maxx-shift(ii);
        newy(idx,ii) = y(1,ii);
    end
else
    newx = x;
    newy = y;
end
% --------------------------------------------------------------------
function varargout = mFileOpenEditor_Callback(h, eventdata, handles, varargin)

num = get(handles.popSel, 'Value');
cur_list = get(handles.list,'String');

for ii =1:length(cur_list)
    cur_commons{ii} = cur_list{ii}(1:2);
    cur_textonly{ii} = cur_list{ii}(3:end);
end
outscr = kv_scripteditor(cur_textonly);
if isempty(outscr), return; end

%%%%%%%% elaborate later @@@@ alsi 14/02/2013
keep_commons = 0;

if ~keep_commons
    for ii = 1:length(outscr)
        list_togo{ii} = ['l-', outscr{ii}];
    end
else
    for ii = 1:length(outscr)
        list_togo{ii} = [cur_commons{ii}, outscr{ii}];
    end
end
set(handles.list, 'String', list_togo, 'Value', 1);
handles.sim{get(handles.popSel, 'Value')}.Script = outscr;
guidata(handles.figure1, handles);
list_Callback(h, [], handles);


% --------------------------------------------------------------------
function varargout = mFileReload_Callback(h, eventdata, handles, varargin)

num = get(handles.popSel, 'Value');
if isfield(handles.sim{num}, 'ScriptName'),
    if ~isempty(handles.sim{num}.ScriptName),
        [name, script] = LoadScript(handles.sim{num}.ScriptName, handles.sim{num}.stype);
        set(handles.list, 'String', script, 'Value', 1);

        guidata(handles.figure1, handles);
        list_Callback(h, [], handles);
    end
end

% --- Executes on button press in pbReCalc.
function pbReCalc_Callback(hObject, eventdata, handles)
ePars_Callback(hObject, [], handles);
% --------------------------------------------------------------------
function varargout = mShsim_Callback(hObject, eventdata, handles, varargin)
mSimShow_Callback(hObject, eventdata, handles);

% --------------------------------------------------------------------
function mCheckOffOn_Callback(hObject, eventdata, handles)
str = {'off', 'on'};
stat = strcmp(get(hObject, 'Checked'), 'on');
set(hObject, 'Checked', str{~stat + 1});
plothis(handles);

function handles = figuresetup(handles)
[fpath,~,e] = fileparts(which('kazan'));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% MENUS
    handles.mFile = uimenu(handles.figure1, 'Tag', '', 'Label', 'File', 'Checked', 'off', 'Separator', 'off', 'Callback', '' );
    handles.AddSim = uimenu(handles.mFile, 'Tag', 'AddSim', 'Label', 'Add sim ...', 'Checked', 'off', 'Separator', 'off', 'Callback', 'simplugin(''AddSim_Callback'',gcbo,[],guidata(gcbo))');
    handles.mFileBaseline = uimenu(handles.mFile, 'Tag', 'mFileBaseline', 'Label', 'Add baseline ...', 'Checked', 'off', 'Separator', 'off', 'Callback', 'simplugin(''AddSim_Callback'',gcbo,[],guidata(gcbo))' );
    handles.mLoad = uimenu(handles.mFile, 'Tag', 'mLoad', 'Label', 'Load ...', 'Checked', 'off', 'Separator', 'off', 'Callback', 'simplugin(''mLoad_Callback'',gcbo,[],guidata(gcbo))' );
    handles.mLoadAll = uimenu(handles.mFile, 'Tag', 'mLoadAll', 'Label', 'Load All...', 'Checked', 'off', 'Separator', 'off', 'Callback', 'simplugin(''mLoad_Callback'',gcbo,[],guidata(gcbo))' );
    handles.mFileCopy = uimenu(handles.mFile, 'Tag', 'mFileCopy', 'Label', 'Copy', 'Checked', 'off', 'Separator', 'off', 'Callback', 'simplugin(''mFileCopy_Callback'',gcbo,[],guidata(gcbo))');
    handles.mFileReload = uimenu(handles.mFile, 'Tag', 'mFileReload', 'Label', 'Reload', 'Checked', 'off', 'Separator', 'off', 'Callback', 'simplugin(''mFileReload_Callback'',gcbo,[],guidata(gcbo))' );
    handles.mSave = uimenu(handles.mFile, 'Tag', 'mSave', 'Label', 'Save ...', 'Checked', 'off', 'Separator', 'off', 'Callback', 'simplugin(''mSave_Callback'',gcbo,[],guidata(gcbo))' );
    handles.mSaveAll = uimenu(handles.mFile, 'Tag', 'mSaveAll', 'Label', 'Save All ...', 'Checked', 'off', 'Separator', 'off', 'Callback', 'simplugin(''mSave_Callback'',gcbo,[],guidata(gcbo))' );
    handles.mFileDelete = uimenu(handles.mFile, 'Tag', 'mFileDelete', 'Label', 'Delete', 'Checked', 'off', 'Separator', 'off', 'Callback', 'simplugin(''butDel_Callback'',gcbo,[],guidata(gcbo))' );
    handles.mDeleteAll = uimenu(handles.mFile, 'Tag', 'mDeleteAll', 'Label', 'Delete All', 'Checked', 'off', 'Separator', 'off', 'Callback', 'simplugin(''butDel_Callback'',gcbo,[],guidata(gcbo))' );
    
    handles.mEdit = uimenu(handles.figure1, 'Tag', 'mEdit', 'Label', 'Edit', 'Checked', 'off', 'Separator', 'off', 'Callback', '' );
    handles.mRefresh = uimenu(handles.mEdit, 'Tag', 'mRefresh', 'Label', 'Refresh', 'Checked', 'off', 'Separator', 'on', 'Callback', 'simplugin(''mRefresh_Callback'',gcbo,[],guidata(gcbo))' );
    handles.mFileChName = uimenu(handles.mEdit, 'Tag', 'mFileChName', 'Label', 'Simulation Name ...', 'Checked', 'off', 'Separator', 'off', 'Callback', 'simplugin(''mFileChName_Callback'',gcbo,[],guidata(gcbo))' );
    handles.mSortPars = uimenu(handles.mEdit, 'Tag', 'mSortPars', 'Label', 'Sort Parameters', 'Checked', 'off', 'Separator', 'on', 'Callback', 'simplugin(''mSortPars_Callback'',gcbo,[],guidata(gcbo))' );
    handles.mMod = uimenu(handles.mEdit, 'Tag', '', 'Label', 'Script Modification', 'Checked', 'off', 'Separator', 'on', 'Callback', '' );
    handles.mFileMoveDown = uimenu(handles.mMod, 'Tag', 'mFileMoveDown', 'Label', 'Add Parameter', 'Checked', 'off', 'Separator', 'on', 'Callback', 'simplugin(''mFileMove_Callback'',gcbo,[],guidata(gcbo))' );
    handles.mFileMoveUp   = uimenu(handles.mMod, 'Tag', 'mFileMoveUp', 'Label', 'Add Parameter', 'Checked', 'off', 'Separator', 'on', 'Callback', 'simplugin(''mFileMove_Callback'',gcbo,[],guidata(gcbo))' );
    handles.mFileDelPar   = uimenu(handles.mMod, 'Tag', 'mFileDelPar', 'Label', 'Add Parameter', 'Checked', 'off', 'Separator', 'on', 'Callback', 'simplugin(''mFileDelPar_Callback'',gcbo,[],guidata(gcbo))' );
    handles.mFileChPar    = uimenu(handles.mMod, 'Tag', 'mFileChPar', 'Label', 'Add Parameter', 'Checked', 'off', 'Separator', 'on', 'Callback', 'simplugin(''mFileChPar_Callback'',gcbo,[],guidata(gcbo))' );
    handles.mFileAddPar   = uimenu(handles.mMod, 'Tag', 'mFileAddPar', 'Label', 'Add Parameter', 'Checked', 'off', 'Separator', 'on', 'Callback', 'simplugin(''mFileAddPar_Callback'',gcbo,[],guidata(gcbo))' );
    
    handles.mFileOpenEditor = uimenu(handles.mEdit, 'Tag', 'mFileOpenEditor', 'Label', 'Open in Editor', 'Checked', 'off', 'Separator', 'off', 'Callback', 'simplugin(''mFileOpenEditor_Callback'',gcbo,[],guidata(gcbo))' );
    handles.mFileCopyClipboard = uimenu(handles.mEdit, 'Tag', 'mFileCopyClipboard', 'Label', 'Copy to Clipboard', 'Checked', 'off', 'Separator', 'off', 'Callback', 'simplugin(''mSave_Callback'',gcbo,[],guidata(gcbo))' );
    
    handles.mOptions = uimenu(handles.figure1, 'Tag', 'mOptions', 'Label', 'Options', 'Checked', 'off', 'Separator', 'off', 'Callback', '' );
    handles.mShsim = uimenu(handles.mOptions, 'Tag', 'mShsim', 'Label', 'Show sim', 'Checked', 'on', 'Separator', 'off', 'Callback', 'simplugin(''mShsim_Callback'',gcbo,[],guidata(gcbo))' );
    handles.mShsumm = uimenu(handles.mOptions, 'Tag', 'mShsumm', 'Label', 'Show sum', 'Checked', 'off', 'Separator', 'off', 'Callback', 'simplugin(''mSimShow_Callback'',gcbo,[],guidata(gcbo))' );
    handles.mShdiff = uimenu(handles.mOptions, 'Tag', 'mShdiff', 'Label', 'Show diff', 'Checked', 'off', 'Separator', 'off', 'Callback', 'simplugin(''mSimShow_Callback'',gcbo,[],guidata(gcbo))' );
    handles.mSimSubBl = uimenu(handles.mOptions, 'Tag', 'mSimSubBl', 'Label', 'Subtract baseline', 'Checked', 'off', 'Separator', 'off', 'Callback', 'simplugin(''mSimShow_Callback'',gcbo,[],guidata(gcbo))' );
    handles.mSimAutoBaseline = uimenu(handles.mOptions, 'Tag', 'mSimAutoBaseline', 'Label', 'Auto baseline', 'Checked', 'off', 'Separator', 'off', 'Callback', 'simplugin(''ZZ_Callback'',gcbo,[],guidata(gcbo))' );
    handles.mColSum = uimenu(handles.mOptions, 'Tag', 'mColSum', 'Label', '"Sum" color', 'Checked', 'off', 'Separator', 'on', 'Callback', 'simplugin(''mSimCol_Callback'',gcbo,[],guidata(gcbo))' );
    handles.mColDiff = uimenu(handles.mOptions, 'Tag', 'mColDiff', 'Label', '"Diff" color', 'Checked', 'off', 'Separator', 'off', 'Callback', 'simplugin(''mSimCol_Callback'',gcbo,[],guidata(gcbo))' );
    handles.mSimChi = uimenu(handles.mOptions, 'Tag', 'mSimChi', 'Label', 'Show Chi^2', 'Checked', 'off', 'Separator', 'off', 'Callback', 'simplugin(''mSimChi_Callback'',gcbo,[],guidata(gcbo))' );
    handles.mSimLarm = uimenu(handles.mOptions, 'Tag', 'mSimLarm', 'Label', 'Shift to Larmor Frequency', 'Checked', 'off', 'Separator', 'off', 'Callback', 'simplugin(''mSimLarm_Callback'',gcbo,[],guidata(gcbo))' );
    handles.mCalc = uimenu(handles.mOptions, 'Tag', 'mCalc', 'Label', 'Calculate', 'Checked', 'off', 'Separator', 'on', 'Callback', 'simplugin(''check_Callback'',gcbo,[],guidata(gcbo))' );
    handles.mSimAgr = uimenu(handles.mOptions, 'Tag', 'mSimAgr', 'Label', 'Apply Source Processing', 'Checked', 'on', 'Separator', 'off', 'Callback', 'simplugin(''mSimAgr_Callback'',gcbo,[],guidata(gcbo))' );
    handles.mOptAddScripts = uimenu(handles.mOptions, 'Tag', 'mOptAddScripts', 'Label', 'Add Scripts to DSC', 'Checked', 'off', 'Separator', 'off', 'Callback', 'simplugin(''mCheckOffOn_Callback'',gcbo,[],guidata(gcbo))' );
    handles.mSimLoadFreq = uimenu(handles.mOptions, 'Tag', 'mSimLoadFreq', 'Label', 'Load Frequency', 'Checked', 'off', 'Separator', 'on', 'Callback', 'simplugin(''mSimLoadFreq_Callback'',gcbo,[],guidata(gcbo))' );
    handles.mSimApplyFreq = uimenu(handles.mOptions, 'Tag', 'mSimApplyFreq', 'Label', 'Auto Load Frequency', 'Checked', 'off', 'Separator', 'off', 'Callback', 'simplugin(''mSimApplyFreq_Callback'',gcbo,[],guidata(gcbo))' );
    handles.mSimLoadRange = uimenu(handles.mOptions, 'Tag', 'mSimLoadRange', 'Label', 'Load Range', 'Checked', 'off', 'Separator', 'off', 'Callback', 'simplugin(''mSimLoadRange_Callback'',gcbo,[],guidata(gcbo))' );
    handles.mIS = uimenu(handles.mOptions, 'Tag', '', 'Label', 'Isotopes (EasySpin)', 'Checked', 'off', 'Separator', 'off', 'Callback', 'isotopes' );
    
    
    handles.Untitled_1 = uimenu(handles.figure1, 'Tag', 'Untitled_1', 'Label', 'About...', 'Checked', 'off', 'Separator', 'off', 'Callback', '' );
    handles.mAboutHelp = uimenu(handles.Untitled_1, 'Tag', 'mAboutHelp', 'Label', 'Help', 'Checked', 'off', 'Separator', 'off', 'Callback', 'simplugin(''mAbout_Callback'',gcbo,[],guidata(gcbo))' );
    handles.mAbout = uimenu(handles.Untitled_1, 'Tag', 'mAbout', 'Label', 'About', 'Checked', 'off', 'Separator', 'off', 'Callback', 'simplugin(''mAbout_Callback'',gcbo,[],guidata(gcbo))' );
    
    %%%%%%%%%%%%%%%%%% CONTEXT MENU
    handles.mContext1 = uicontextmenu('Parent', handles.figure1);%(handles.figure1, 'Tag', 'mContext1'); %, 'Callback', 'simplugin(''mContext1_Callback'',gcbo,[],guidata(gcbo))' );
    handles.mCont1AddRange = uimenu(handles.mContext1, 'Tag', 'mCont1AddRange', 'Label', 'Add Range', 'Checked', 'off', 'Separator', 'off', 'Callback', 'simplugin(''mCont1AddRange_Callback'',gcbo,[],guidata(gcbo))' );
    
    handles.mContext = uicontextmenu('Parent', handles.figure1);%(handles.figure1, 'Tag', 'mContext', 'Callback', '' );
    %       handles.mContOpt = uicontextmenu(handles.mContext, 'Tag', 'mContOpt', 'Label', 'Optimise', 'Checked', 'off', 'Separator', 'off', 'Callback', 'simplugin(''mContOptSet_Callback'',gcbo,[],guidata(gcbo))' );
    %       handles.mContNotOpt = uicontextmenu(handles.mContext, 'Tag', 'mContNotOpt', 'Label', 'Do not optimise', 'Checked', 'off', 'Separator', 'off', 'Callback', 'simplugin(''mContOptSet_Callback'',gcbo,[],guidata(gcbo))' );
    handles.mMoveUp = uimenu(handles.mContext, 'Tag', 'mMoveUp', 'Label', 'MoveUp', 'Checked', 'off', 'Separator', 'on', 'Callback', 'simplugin(''mFileMove_Callback'',gcbo,[],guidata(gcbo))' );
    handles.mMoveDown = uimenu(handles.mContext, 'Tag', 'mMoveDown', 'Label', 'Move Down', 'Checked', 'off', 'Separator', 'off', 'Callback', 'simplugin(''mFileMove_Callback'',gcbo,[],guidata(gcbo))' );
    handles.mConMakeCom = uimenu(handles.mContext, 'Tag', 'mConMakeCom', 'Label', 'Make it common', 'Checked', 'off', 'Separator', 'on', 'Callback', 'simplugin(''mConComm_Callback'',gcbo,[],guidata(gcbo))' );
    handles.mConMakeLoc = uimenu(handles.mContext, 'Tag', 'mConMakeLoc', 'Label', 'Make it local', 'Checked', 'off', 'Separator', 'off', 'Callback', 'simplugin(''mConComm_Callback'',gcbo,[],guidata(gcbo))' );
    handles.mContAddPar = uimenu(handles.mContext, 'Tag', 'mContAddPar', 'Label', 'Add Parameters', 'Checked', 'off', 'Separator', 'on', 'Callback', 'simplugin(''mFileAddPar_Callback'',gcbo,[],guidata(gcbo))' );
    handles.mContChPar = uimenu(handles.mContext, 'Tag', 'mContChPar', 'Label', 'Delete Parameter', 'Checked', 'off', 'Separator', 'off', 'Callback', 'simplugin(''mFileDelPar_Callback'',gcbo,[],guidata(gcbo))' );
    handles.mContDelPar = uimenu(handles.mContext, 'Tag', 'mContDelPar', 'Label', 'Change Parameter', 'Checked', 'off', 'Separator', 'off', 'Callback', 'simplugin(''mFileChPar_Callback'',gcbo,[],guidata(gcbo))' );
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SIMULATION
    handles.frame1 = uicontrol(handles.figure1, 'Style', 'frame', 'Tag', 'frame1', 'Units', 'pixels');
    handles.popSel = uicontrol(handles.figure1, 'Style', 'popupmenu', 'Tag', 'popSel', 'String', '1', 'Value', 1 , 'Units', 'pixels', 'Callback', 'simplugin(''popSel_Callback'',gcbo,[],guidata(gcbo))');
    handles.butCol = uicontrol(handles.figure1, 'Style', 'pushbutton', 'Tag', 'butCol',...
        'String', '','Units', 'pixels', 'Callback', 'simplugin(''mSimCol_Callback'',gcbo,[],guidata(gcbo))');
    handles.list = uicontrol(handles.figure1, 'Style', 'listbox', 'Tag', 'list', 'uicontextmenu', handles.mContext, ...
        'Units', 'pixels', 'Callback', 'simplugin(''list_Callback'',gcbo,[],guidata(gcbo))', 'BackgroundColor', [1, 1, 1]);
    
    set(handles.figure1, 'uicontextmenu', handles.mContext);
    set(handles.list, 'String', {'l-Sys.S = 0.5', ...
        'l-Sys.g(1) = 2.1', ...
        'l-Sys.g(2) = 2', ...
        'l-Sys.g(3) = 1.98', ...
        'l-Sys.lw = 3', ...
        'l-Exp.Range = [330,370]', ...
        'l-Sys.gStrain(1) = 0', ...
        'l-Sys.gStrain(2) = 0', ...
        'l-Sys.gStrain(3) = 0', ...
        'l-Exp.mwFreq = 9.77', ...
        'l-Exp.Harmonic = 1', ...
        'l-Opt.Sim = ''simpleEPR''', ...
        'l-Shift.x = 0', ...
        'l-Shift.s = 0', ...
        'l-Shift.Scale = 1'}, ...
        'Value', 1);
    
    handles.ePars = spinbox(handles.figure1, 'Tag', 'ePars', 'Value', 0.5, 'Step', 0.1,'Units', 'pixels', 'Callback', 'simplugin(''ePars_Callback'',[],[], guidata(gcf))');
    handles.popup = uicontrol(handles.figure1, 'Style', 'popupmenu', 'Tag', 'popup',...
        'String', { '1e-6', '1e-5', '1e-4', '1e-3', '1e-2', '0.1', '1', '10', '100', '1e+3', '1e+4', '1e+5', '1e+6'}, ...
        'Value', 6, 'Units', 'pixels', 'Callback', 'simplugin(''popup_Callback'',gcbo,[],guidata(gcbo))');
    handles.text4 = uicontrol(handles.figure1, 'Style', 'text', 'Tag', 'text4', 'String', 'Simulation', 'Units', 'pixels');
    
    %%%%%%%%%%%%%%%%%%% "FIT"
    handles.frame2 = uicontrol(handles.figure1, 'Style', 'frame', 'Tag', 'frame2', 'Units', 'pixels', 'Callback', '');
    
    
    handles.pmOptMethod = uicontrol(handles.figure1, 'Style', 'popupmenu', 'Tag', 'pmOptMethod', 'Units', 'pixels', 'Callback', 'simplugin(''ZZ_Callback'',gcbo,[],guidata(gcbo))');
    set(handles.pmOptMethod, 'String', {'None', ...
        'Fit Peak-To-Peak', ...
        'Fit Amplitude', ...
        'Fit Integral', ...
        'Fit Double Integral', ...
        'Fix Double Integral'}, ...
        'Value', 1);
    %         'Default optimization', ...
    %         '10 Iterations', ...
    %         '100 Iterations', ...
    %         '1.0E-6 tolerance', ...
    %         '1.0E-8 tolerance'}', ...
    
    handles.lbOpt = uicontrol(handles.figure1, 'Style', 'listbox', 'Tag', 'lbOpt', 'String', { 'Opm.Amp = 1', 'Opm.wf = 0', 'Opm.dys = 0'}, 'Value', 1, 'Units', 'pixels', 'Callback', 'simplugin(''lbOpt_Callback'',gcbo,[],guidata(gcbo))');
    handles.pmOpt = uicontrol(handles.figure1, 'Style', 'popupmenu', 'Tag', 'pmOpt', 'String', { '1e-6', '1e-5', '1e-4', '1e-3', '1e-2', '0.1', '1', '10', '100', '1e+3', '1e+4', '1e+5', '1e+6', '1e+7', '1e+8'},'Value', 6, 'Units', 'pixels', 'Callback', 'simplugin(''pmOpt_Callback'',gcbo,[],guidata(gcbo))');
    handles.eOpt = spinbox(handles.figure1, 'Tag', 'eOpt', 'Value', 1, 'Step', 0.1,'Units', 'pixels', 'Callback', 'simplugin(''eOpt_Callback'',[],[], guidata(gcf))');
    handles.text3 = uicontrol(handles.figure1, 'Style', 'text', 'Tag', 'text3', 'String', 'fit', 'Units', 'pixels');
    
    handles.pbReCalc = uicontrol(handles.figure1, 'Style', 'pushbutton', 'Tag', 'pbReCalc', 'String', '',  'Units', 'pixels', 'Callback', 'simplugin(''pbReCalc_Callback'',gcbo,[],guidata(gcbo))');
    handles.pbCalcAll = uicontrol(handles.figure1, 'Style', 'pushbutton', 'Tag', 'pbCalcAll', 'String', '', 'Units', 'pixels', 'Callback', 'simplugin(''calc_all'', guidata(gcbo))');
    handles.butCalc = uicontrol(handles.figure1, 'Style', 'togglebutton', 'Tag', 'butCalc', 'String', '', 'Units', 'pixels', 'Callback', 'simplugin(''butCalc_Callback'',gcbo,[],guidata(gcbo))');
    % 	handles.slOpt = uicontrol(handles.figure1, 'Style', 'slider', 'Tag', 'slOpt', 'String', '{, ''}', 'Value', ' ', 'Units', 'pixels', 'Callback', 'simplugin('slOpt_Callback',gcbo,[],guidata(gcbo))');
    % 	handles.slPars = uicontrol(handles.figure1, 'Style', 'slider', 'Tag', 'slPars', 'String', '{, ''}', 'Value', ' ', 'Units', 'pixels', 'Callback', 'simplugin('slPars_Callback',gcbo,[],guidata(gcbo))');
    
    handles.tbMain = uitoolbar(handles.figure1);
    load([fpath, filesep, 'ico.mat']);
    handles.ptChange = uipushtool(handles.tbMain, 'TooltipString', 'Change script by loading file', 'Tag', 'ptChange', ...
        'CData', ico.change, 'ClickedCallback', 'simplugin(''mLoad_Callback'',gcbo,[],guidata(gcbo))');
    handles.ptOpen = uipushtool(handles.tbMain, 'TooltipString', 'Add script from file', 'Tag', 'ptOpen', ...
        'CData', ico.load, 'ClickedCallback', 'simplugin(''AddSim_Callback'',gcbo,[],guidata(gcbo))');
    
    handles.ptSave = uipushtool(handles.tbMain, 'TooltipString', 'Save current script', 'Tag', 'ptSave',...
        'CData', ico.save, 'ClickedCallback', 'simplugin(''mSave_Callback'',gcbo,[],guidata(gcbo))');
    
    handles.ptCopy = uipushtool(handles.tbMain, 'TooltipString', 'Duplicate Script', 'Tag', 'ptCopy', ...
        'CData', ico.copy, 'ClickedCallback', 'simplugin(''mFileCopy_Callback'',gcbo,[],guidata(gcbo))');
    
    handles.ptEdit = uipushtool(handles.tbMain, 'TooltipString', 'Change Script in editor', 'Tag', 'ptEdit', ...
        'CData', ico.edit, 'ClickedCallback', 'simplugin(''mFileOpenEditor_Callback'',gcbo,[],guidata(gcbo))');
    
    handles.ptReload = uipushtool(handles.tbMain, 'TooltipString', 'Reload', 'Tag', 'ptReload', ...
        'CData', ico.reload, 'ClickedCallback', 'simplugin(''mFileReload_Callback'',gcbo,[],guidata(gcbo))');
    
    handles.ptClipboard = uipushtool(handles.tbMain, 'TooltipString', 'Copy to Clipboard', 'Tag', 'ptClipboard',...
        'CData', ico.clipboard,'ClickedCallback', 'simplugin(''mSave_Callback'',gcbo,[],guidata(gcbo))');
    
    handles.mIsotopes = uimenu(handles.mOptions, 'Label', 'Isotopes (EasySpin)', 'Callback', 'isotopes');
    set(handles.butCalc,  'CData', ico.calcerr, 'BackgroundColor', [0, 0, 0]);
    set(handles.pbReCalc,  'CData', ico.recalculate, 'BackgroundColor', [0, 0, 0]);
    set(handles.pbCalcAll,  'CData', ico.calcall, 'BackgroundColor', [0, 0, 0]);
    handles.hh = 0;
    handles.ploter = 0;
    handles.num = 1;
    handles.div = 0;
    handles.Spec = {};
    handles.sim{1}.Script = get(handles.list, 'String');
    handles.sim{1}.x = [];
    handles.sim{1}.y = [];
    handles.sim{1}.stype = 'sim';
    handles.Color =   [1, 0, 0; 1, .6, 0; 0, 1, 0;...
        0, 0, 0; 1, 0, 1; 0, 1, 1; 1, 1, 0];
    handles.SumCol = [1, 0, 1];
    handles.DiffCol = [0, 1, 0];
    handles.CalcType = 'Pepper';
    set(handles.mFileMoveUp, 'Accelerator', 'P');
    set(handles.mFileMoveDown, 'Accelerator', 'O');
    set(handles.mLoad, 'Accelerator', 'L');
    set(handles.mSave, 'Accelerator', 'S');
    set(handles.mCalc, 'Accelerator', 'C');
    set(handles.AddSim, 'Accelerator', 'A');
    set(handles.mFileBaseline, 'Accelerator', 'B');
    set(handles.mSimChi, 'Accelerator', 'H');
    set(handles.butCol, 'BackgroundColor', handles.Color(1, :));
    set(handles.ePars, 'UserData', ' ');
    set(handles.mSimAgr, 'Checked', 'on');
    handles.process = 1;
    handles.applyfreq = 0;
    handles.Larmfreq  = 0;
    handles.DataRange.bracket = [];
    handles.DataRange.idx = [];
