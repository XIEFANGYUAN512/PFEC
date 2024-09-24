function [Y_pred, obj] = PFEC(alpha,A0,baseY,c)
% The main function of SWEC
%B = Y*theta*Y'
%
[~,M] = size(baseY);
lambda = 5;
cm = zeros(M,1);
for i = 1:M
    [~,cm(i)] = size(baseY{i});
    temp1 = rand(cm(i),1);
    temp1 = temp1/sum(temp1);
    theta{i} = diag(temp1);
    B{i} = baseY{i}*theta{i}*baseY{i}';
    P{i} = (baseY{i}'*baseY{i}).*(baseY{i}'*baseY{i});
end
for iter = 1:30
    %update A
    [Y_pred, A] = CLR_updateA(alpha, B, c, lambda, A0);
    
    %update theta
    for i = 1:M
        Q{i} = diag(baseY{i}'*A*baseY{i}); 
        [theta_diag, ~]=SimplexQP_acc(P{i}, Q{i}, diag(theta{i}));
        theta{i} = diag(theta_diag);
    end
    
    %update alpha
    for i = 1:M
        B{i} = baseY{i}*theta{i}*baseY{i}';
    end
    for i = 1:M
        alpha(i) = 0.5/(norm(A-B{i},'fro'));
    end
    
    %calculate the objective function value
    obj_t = 0;
    for i = 1:M
        obj_t = obj_t + norm(A-B{i},'fro');
    end
    obj(iter) = obj_t;
end


end

