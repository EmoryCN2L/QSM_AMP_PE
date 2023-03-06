function opt = setoptdefaults(opt,o)
%%this function set the default values of options stored in a cell
%function opt = setoptdefaults(opt,o)
%opt, input option structure
%o, a cell containing pairs of option name and default values
%opt, output
%Example:
%opt = setoptdefaults(opt,{'B0',3,...
%                    'TE',30e-3,...
%                    'phasescale',1,...
%                    'tols',1e-t,...
%                    'filterlength',6,...
%                    'coarselevel',3});
%Deqiang Qiu. qiudeqiang@gmail.com
%Stanford University. 2011

if(isempty(opt))
    opt = struct;
end

for i = 1:(numel(o)/2)
    if(~isfield(opt,o{i*2-1}))
        opt.(o{i*2-1}) = o{i*2};
    end
end
return


