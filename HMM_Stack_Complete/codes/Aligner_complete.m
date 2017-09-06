function [core_return,figs] = Aligner_complete(files,data,rhos,core_param,stack_param)

% Criteria of convergence
iterTol = 0;
iterMax = 10; 

% Initial Parameter Values from sedemenation rate
done = false;

P = length(files{1});
PMatrix_Total = cell(P,1);

iter = 0; 
while ~done
    iter = iter + 1;
    old_core_param = core_param;
    new_core_param = core_param;
    
    parfor index = 1:P
        % Forward Algorithm
        [fMatrix,~] = forward_algorithm_complete(data,stack_param,core_param,index,rhos);
        tt = ['Forward algorithim for the core ',data(index).name,' in iteration ',num2str(iter),' is done.'];
        disp(tt);
        
        % Backward Sampling Algorithm
        [BMatrix,max_T] = backward_algorithm_complete(data,stack_param,core_param,index,rhos);
        tt = ['Backward algorithim for the core ',data(index).name,' in iteration ',num2str(iter),' is done.'];
        disp(tt);
        
        % Computing Posterior:
        [PMatrix,QMatrix] = posterior_complete(fMatrix,BMatrix,max_T);
        tt = ['Computing the posterior for the core ',data(index).name,' in iteration ',num2str(iter),' is done.'];
        disp(tt);
        
        % Updating Core Parameters:
        new_core_param(index) = core_parameter_update(fMatrix,BMatrix,data,index,PMatrix,QMatrix,core_param,stack_param);
        tt = ['Updating the core parameters for the core ',data(index).name,' in iteration ',num2str(iter),' is done.'];
        disp(tt);

        PMatrix_Total{index} = PMatrix';
        
    end
    
    core_param = new_core_param;
    
    % Termination Criterion:
    Diff = Diff_SQR(old_core_param,core_param);
    if (iter >= iterMax) || (Diff < iterTol)
        done = true;
    end 
end

core_return = struct('name',cell(P,1),'Posterior',cell(P,1),'shift',cell(P,1),'R',cell(P,1),'alpha',cell(P,1),'beta',cell(P,1),'phi',cell(P,1),'psi',cell(P,1),'eta',cell(P,1),'epsilon',cell(P,1));
for p = 1:P
    core_return(p).name = core_param(p).name;
    core_return(p).shift = core_param(p).shift;
    core_return(p).R = core_param(p).R;
    core_return(p).alpha = core_param(p).alpha;
    core_return(p).beta = core_param(p).beta;
    core_return(p).phi = core_param(p).phi;
    core_return(p).psi = core_param(p).psi;
    core_return(p).eta = core_param(p).eta;
    core_return(p).epsilon = core_param(p).epsilon;
    core_return(p).Posterior = PMatrix_Total{p};
end

figs = plotting_results_complete(PMatrix_Total,data,stack_param,files);


cd ..

end
