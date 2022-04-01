function spinboxcb(hands, obj)

if ~isa(hands, 'spinbox'), return; end

param = get(hands.push2, 'userdata');
if strcmp(param.Enable, 'off'), return; end
donothing = 0;
switch obj
    case 'edit'
            param.Value = str2num(get(hands.edit, 'String'));
            param.String = get(hands.edit, 'String');
            if ~isempty(param.Value)
                set(hands, 'EnableButtons', 'on');
            else
                
                set(hands, 'EnableButtons', 'off');
            end
            sss = builtin('get', hands.edit, 'String');
            
    case 'push2' % down
        if ~ischar(param.Value)
            num = param.Value + param.Step;
            param.Value = num;
            param.String = num2str(param.Value);
        end
        builtin('set', hands.edit, 'String', num2str(num));
    case 'push1' % up
        if ~ischar(param.Value)
            num = param.Value - param.Step;
            param.Value = num;
            param.String = num2str(param.Value);
        end
        builtin('set', hands.edit, 'String', num2str(num));
    case 'keypress'
        donothing = 1;
        key = double(builtin('get', gcf, 'CurrentCharacter'));
        if ~isempty(key)
            if key==30 %|| key==29 % up key or right
                if ~ischar(param.Value)
                    num = param.Value + param.Step;
                    param.Value = num;
                    param.String = num2str(num);
                end
                builtin('set', hands.edit, 'String', num2str(num));
                donothing = 0;
            elseif key==31 %|| key==28  % down key or left
                if ~ischar(param.Value)
                    num = param.Value - param.Step;
                    param.Value = num;
                    param.String = num2str(num);
                end
                builtin('set', hands.edit, 'String', num2str(num));
                donothing = 0;
            end
       
        end
end
if ~donothing  
builtin('set', hands.push2, 'userdata', param);
  
    if ~isempty(param.Callback)
        %     global ed;
        ed = hands;
        eval(param.Callback)
    end
end

function hhh = gcbo
% just for case, when in callback is 'gcbo' function
% global ed;
hhh = [];%builtin('gcbo');

function hhh = gco
% just for case, when in callback is 'gco' function
% global ed;
hhh = builtin('gco');