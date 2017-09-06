stack_name = 'Prob_Stack_short';

% If you just want to align a core to the stack with both dO18 and radiocarbon:
% [core_return,figs] = Aligner('record_summary.txt',stack_name,'dO18','radiocarbon');
% If you have dO18 data only:
[core_return,figs] = Aligner('record_summary.txt',stack_name,'dO18');

% TODO: If you want to learn the stack:
% [stack_return,figs] = Learner('record_summary.txt','dO18','radiocarbon');