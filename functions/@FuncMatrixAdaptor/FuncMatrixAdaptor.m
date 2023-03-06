function obj = FuncMatrixAdaptor(mtimes,t_mtimes)
%This class adapts function handles to matrix product operation through
%function override
%
%function obj = func2MatrixAdpator(mtimes,t_mtimes)
%Deqiang Qiu
%
obj.adjoin = 0;

obj.Ops.mtimes = mtimes;
obj.Ops.t_mtimes = t_mtimes;

obj = class(obj,'FuncMatrixAdaptor');

end


