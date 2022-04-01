function fittingbox_func(varargin)
% fittingbox_func Application M-file ONLY for fittingbox.m
% Construct GUI for selecting functions
%
% input argument: handles of fittingbox uicontrols
%
% Last Modified 23.02.2004
% Silakov Alexey 23-feb-2004, homework

if nargin == 0
elseif nargin == 1, % launch dialog
    handles = varargin{1};
    fgFuncList = figure('MenuBar', 'none', 'Resize', 'on');
    set(fgFuncList, 'Tag', 'FuncList');
    set(fgFuncList, 'Name', 'Select function');
    set(fgFuncList, 'PaperUnits', 'normalized');
    set(fgFuncList, 'Position', [103 200 500  300]);
    set(fgFuncList, 'UserData', handles);
    set(fgFuncList, 'CloseRequestFcn', 'fittingbox_func(''Close_Callback'', gcbo)');
    
    stPar = uicontrol('Parent', fgFuncList, 'Tag', 'stPar');
    set(stPar, 'Style', 'text');
    set(stPar, 'Units', 'normalized');
    set(stPar, 'Position', [0.045 0.15 0.4 0.3]);
    set(stPar, 'String', 'Parameters');
    
    lbFunc = uicontrol('Parent', fgFuncList, 'Tag', 'lbFunc');
    set(lbFunc, 'Style', 'listbox');
    set(lbFunc, 'Units', 'normalized');
    set(lbFunc, 'Position', [0.01 0.47 0.45 0.48]);
    set(lbFunc, 'Callback', 'fittingbox_func(''lbFunc_Callback'', gcbo)');
    functions = {'polyfit(x, y, n)';...
            'exppolyfit(x, y, n)';...
            'kv_expdata(t1, t2, x0, x, y)'; ...
            'a*exp(-(x-xo)/t) + yo';...
            'a1*exp(-(x-xo)/t1) + a2*exp(-(x-xo)/t2) + yo';...
            'a*(1-exp(-(x-xo)/t)) + yo';...
            'a1*(1-exp(-(x-xo)/t1)) + a2*(1-exp(-(x-xo)/t2)) + yo';...
            'a*sin((x - xo)*b) + yo';...
            'a*cos((x - xo)*b) + yo';...
            'a*lshape(x,xo,fwhh,diff,alpha)';...
            'a + b*x + c*x.^2'; ...
            'fit_deer_gaus(x, r0, dr, shi, amp)'};
    parameters.par = {{'n = 3'};...
            {'n = 3'};...
            {'t1 = 1'; 't2 = 5'; 'x0 = 0'}; ...
            {'a = 1'; 'xo = -0.5'; 't = 2'; 'yo = 0.5'};...
            {'a1 = -1'; 'xo = -0.5'; 't1 = 1'; 'yo = 0';...
                'a2 = 1'; 't2 = 2'};...
            {'a = 1'; 'xo = -0.5'; 't = 2'; 'yo = 0.5'};...
            {'a1 = -1'; 'xo = -0.5'; 't1 = 1'; 'yo = 0';...
                'a2 = 1'; 't2 = 2'};...
            {'a = 1'; 'xo = 0'; 'b = 3'; 'yo = 0.5'};...
            {'a = 1'; 'xo = 0'; 'b = 3'; 'yo = 0.5'};...
            {'a = 1'; 'xo = 5'; 'fwhh = 2'; 'diff = 1'; 'alpha = 1'};...
            {'a = 1'; 'b = -10'; 'c = 3'};...
            {'r0 = 4.3'; 'dr = 0.3'; 'shi = 0'; 'amp = 1'}};
    
    parameters.addparam(1) = {'x = ([1:100]/10); y = (1 - 0.2*rand(1, 100)/2).*(1 - x + 0.5*x.^2);'};
    parameters.addparam(2) = {'x = ([1:100]/10); y = (1 - 0.2*rand(1, 100)/2).*(2*exp(-x) + 5*exp(-x/5) + 1);'};
    parameters.addparam(3) = {'x = ([1:100]/10)''; y = (1 - 0.2*rand(100, 1)/2).*(2*exp(-x) + 5*exp(-x/5));'};
    parameters.addparam(4:11) = {'x = ([1:100]/10)'';'};
    parameters.addparam(12) = {'x = 0:10:2000; y = fit_deer_gaus(x, 4.3, 0.3, 0, 1, 1)+0.002*(rand(201, 1)-0.5);'};
    
    parameters.comment(1) = {'no fitting parameters'};
    set(lbFunc, 'String', functions);
    set(lbFunc, 'UserData', parameters);
    
    
    pbOK = uicontrol('Parent', fgFuncList, 'Tag', 'pbOK');
    set(pbOK, 'Style', 'pushbutton');
    set(pbOK, 'Units', 'normalized');
    set(pbOK, 'Position', [0.01 0.02 0.2 0.1]);
    set(pbOK, 'Callback', 'fittingbox_func(''pbOK_Callback'', gcbo)');
    set(pbOK, 'String', 'OK');
    
    pbCancel = uicontrol('Parent', fgFuncList, 'Tag', 'pbCancel');
    set(pbCancel, 'Style', 'pushbutton');
    set(pbCancel, 'Units', 'normalized');
    set(pbCancel, 'Position', [0.25 0.02 0.2 0.1]);
    set(pbCancel, 'Callback', 'fittingbox_func(''pbCancel_Callback'', gcbo)');
    set(pbCancel, 'String', 'Cancel');
    
    axFunc = axes('Parent', fgFuncList); 
    set(axFunc, 'Tag', 'axFunc');
    set(axFunc, 'Units', 'normalized');
    set(axFunc, 'Position', [0.53 0.15 0.45 0.80]);
    
elseif ischar(varargin{1}) % INVOKE NAMED SUBFUNCTION OR CALLBACK
    
    try
        if (nargout)
            [varargout{1:nargout}] = feval(varargin{:}); % FEVAL switchyard
        else
            feval(varargin{:}); % FEVAL switchyard
        end
    catch
        disp(lasterr);
    end
end

%---------------------------------------------------
function pbOK_Callback(h)
fighand = get(h, 'Parent');
hhh = get(fighand, 'UserData');
hand = guidata(hhh.figure1);

listhand = findobj('Tag','lbFunc');
str = get(listhand, 'String');
num = get(listhand, 'Value');
par = get(listhand, 'UserData');
param = par.par{num};
restr{1} = strcat('f = ', str{num});
for count = 1:size(param)
    restr{end + 1} = param{count};
end

set(hand.list, 'UserData', restr);
set(fighand, 'UserData', str{num});
uiresume(hand.figure1);
delete(fighand);

%---------------------------------------------------
function pbCancel_Callback(h)
fighand = get(h, 'Parent');
hhh = get(fighand, 'UserData');
hand = guidata(hhh.figure1);
set(fighand, 'UserData', []);
uiresume(hand.figure1);
delete(fighand);

%---------------------------------------------------
function lbFunc_Callback(h)
str = get(h, 'String');
num = get(h, 'Value');
par = get(h, 'UserData');

func = str{num};
param = par.par{num};
addp = par.addparam{num};
 
eval(addp);
for count = 1:size(param, 1)
    eval([param{count}, ';']);
end
if num == 1
    f(:, 1) = polyval(polyfit(x', y', n), x');
    f(:, 2) = y';
elseif num == 2
 
    pl = polyfit(x', log(y'), n);
    f(:, 1) = exp(polyval(pl, x'));
    f(:, 2) = y';
else
    eval(['f =', func, ';']);
end
if sum(num == [3, 12]), f(:, 1) = f; f(:, 2) = y; end

hands = get(get(h, 'Parent'), 'Children');
axhand = findobj(hands, 'Type', 'axes');
txthand = findobj(hands, 'Tag', 'stPar');
set(txthand, 'String', param);
col = [1, 0, 0; 0, 0, 1];
for i = 1:size(f, 2)
    plot(x(:), f(:, i), 'Color', col(i, :), 'Parent', axhand); hold on;
end
hold off;
axis auto;

%---------------------------------------------------
function Close_Callback(h)
try
    hand = guidata(get(h, 'UserData'));
    uiresume(hand.figure1);
end
delete(h);