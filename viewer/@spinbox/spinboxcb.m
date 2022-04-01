function spinboxcb(hands, obj)

if ~isa(hands, 'spinbox'), return; end

param = get(hands.push2, 'userdata');
if strcmp(param.Enable, 'off'), return; end
switch obj
    case 'edit'
            param.Value = str2num(get(hands.edit, 'String'));
            if ~isempty(param.Value)
                set(hands, 'EnableButtons', 'on');
            else
                param.Value = get(hands.edit, 'String');
                set(hands, 'EnableButtons', 'off');
            end
            sss = builtin('get', hands.edit, 'String')
    case 'push2' % down
        if ~ischar(param.Value)
            num = param.Value + param.Step;
            param.Value = num;
        end
        builtin('set', hands.edit, 'String', num2str(num));
    case 'push1' % up
        if ~ischar(param.Value)
            num = param.Value - param.Step;
            param.Value = num;
        end
        builtin('set', hands.edit, 'String', num2str(num));
end
builtin('set', hands.push2, 'userdata', param);


if ~isempty(param.Callback)
    global ed;
    ed = hands;
    eval(param.Callback)
end

function hhh = gcbo
% just for case, when in callback is 'gcbo' function
global ed;
hhh = ed;

function hhh = gco
% just for case, when in callback is 'gco' function
global ed;
hhh = ed;