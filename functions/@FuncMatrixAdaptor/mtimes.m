function res = mtimes(obj,x)

if(obj.adjoin==0)
    if(isempty(obj.Ops.mtimes))
        error('obj.Ops.mtimes is empty');
    else
        res = obj.Ops.mtimes(x);
        return;
    end
else
    if(isempty(obj.Ops.t_mtimes))
        error('obj.Ops.t_mtimes is empty')
    else
        res = obj.Ops.t_mtimes(x);
        return;
    end
end

return

        