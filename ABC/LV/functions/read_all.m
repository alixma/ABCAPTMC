function [output_out, n_out, TM_out, stats_out, rates_out] = read_all(dirname, name, ntimes, output_in, n_in, TM_in, stats_in, rates_in)
% read everything and add to existing results

[output_out, n_out, TM_out, stats_out, rates_out] = read_one(dirname, name, ntimes);
% in case we already have some output and want to add a previous run
if nargin>3 
    output_out = [output_in output_out]; TM_out = [TM_in TM_out];
    n_out = [n_in n_out]; stats_out = cat(3, stats_in, stats_out);
    rates_out = cat(3, rates_in, rates_out);
end




end