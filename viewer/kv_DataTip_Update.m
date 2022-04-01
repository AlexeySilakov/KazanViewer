% --------------------------------------------------------------------
%for ML>7. Set to make only 
function [txt] = kv_DataTip_Update(obj,event_obj)
    % Display 'Time' and 'Amplitude'
pos = get(event_obj,'Position');
txt = {num2str(pos(1))};
