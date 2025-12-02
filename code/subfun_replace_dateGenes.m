function [ OutGeneSym ] = subfun_replace_dateGenes( DateSym )

tmploc = find(DateSym=='-');
if isempty(tmploc)
    tmploc = find(DateSym=='.');
end
switch DateSym(tmploc+1:end)
    case 'Mar'
        OutGeneSym = ['MARCH' DateSym(1:tmploc-1)];
    case 'Jun'
        OutGeneSym = ['JUN' DateSym(1:tmploc-1)];
    case 'Sep'
        OutGeneSym = ['SEPT' DateSym(1:tmploc-1)];
    case 'Dec'
        OutGeneSym = ['DECR' DateSym(1:tmploc-1)];
end

end

