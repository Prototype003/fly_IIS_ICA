function [ sigs, p_thresh ] = fdr_correct(ps, q)
% Applies FDR correction to a vector of p values

[p_thresh, p_thresh_nonpara] = FDR(ps, q);

if isempty(p_thresh) == 0
    sigs = ps < p_thresh;
else
    sigs = zeros(length(ps), 1);
end

end

