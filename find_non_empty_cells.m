function [inds,inds2]=find_non_empty_cells(x)
if nargout<2
    inds=find(~cellfun('isempty', x));
else
    [inds,inds2]=find(~cellfun('isempty', x));
end