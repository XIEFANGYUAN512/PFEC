% The implement of "Parameter-free ensemble clustering with dynamic
% weighting mechanism". If it is useful, please cite:
% @article{XIE2024110389,
% title = {Parameter-free ensemble clustering with dynamic weighting mechanism},
% journal = {Pattern Recognition},
% volume = {151},
% pages = {110389},
% year = {2024},
% issn = {0031-3203},
% doi = {https://doi.org/10.1016/j.patcog.2024.110389},
% url = {https://www.sciencedirect.com/science/article/pii/S0031320324001407},
% author = {Fangyuan Xie and Feiping Nie and Weizhong Yu and Xuelong Li},
% }

clear;close all;
addpath(genpath(pwd));

Dataset = load('bp_wine_unikmcsqrtn.mat');
Y = Dataset.Y;baseY_o = Dataset.baseCls;
[n,~] = size(Y);c = length(unique(Y));repeat_max = 10;
%% ensemble clusters generation
for repeat_num = 1:repeat_max

    [~,M_all] = size(baseY_o);M = 10;
    baseY = cell(1,M);cm = zeros(M,1);
    temp = randperm(M_all);
    temp = temp(1:M);
    A0 = zeros(n);

    alpha = rand(1,M);
    alpha = alpha/sum(alpha);
    for i = 1:M
        baseY_nM(:,i) = baseY_o(:,temp(i));
        baseY{i} = n2nc(baseY_nM(:,i));
        A0 = A0 + baseY{i}*(baseY{i}'*baseY{i})^(-1)*baseY{i}';%initialize A
    end
    A0 = A0/M;

    % the performance of base clusterings
    for num = 1:M
        base_acc(num,:) = ClusteringMeasure_All(Y,nc2n(baseY{num}));
    end
    base_accs{repeat_num} = base_acc;

    %% The PFEC
    [Y_pfec, obj] = PFEC(alpha,A0,baseY,c);
%             figure();plot(obj);
    PFECs.result(repeat_num,:) = ClusteringMeasure_All(Y,Y_pfec);
end
mean(PFECs.result)
