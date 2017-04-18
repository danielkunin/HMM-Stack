% [Confidence_Band,core_param,stack_param,samples,sites,data] = Aligner_del_O18('record_summary.txt');
% [Confidence_Band,core_param,stack_param,samples,sites,data] = Aligner_Intervals('record_summary.txt');
% [Confidence_Band,core_param,stack_param,samples,sites,data] = Aligner_Pseudodata('record_summary.txt');
[Confidence_Band,core_param,stack_param,samples,sites,data] = Aligner_Complete('record_summary.txt');


% plotting_results_del_O18;
% plotting_results_intervals;
% plotting_results_pseudodata;
plotting_results_complete;