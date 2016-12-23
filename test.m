max_iters = 500; conv_tol = 1e-6; total_iters = 100;
min_dim = 2; max_dim = 5;
range_dim = max_dim-min_dim;
n = 10;
prior = [-1 1 10];
valid_undirected = [];
I_undirected = [];
ifconverged_undirected = [];
valid_directed = [];
I_directed = [];
ifconverged_directed = [];
Stats_name = {'# of edges','# of connected components',...
    'diameter(the longest shortest path)','the average shortest path','radius(minimum vertex eccentricity)',...
    '# of size 3 loop','# of size 4 loop','# of independent loop','average factor node dimension(average cardinality)',...
    'maximum factor node dimension(average cardinality)',};
Stats_undirected = zeros(length(Stats_name),total_iters);
Stats_directed = zeros(length(Stats_name),total_iters);
for i=1:total_iters
    fprintf('%d iteration\n', i);
    facG_root_undirected.var=[];
    facG_root_undirected.fac=[];
    facG_root_directed.var=[];
    facG_root_directed.fac=[];
    while(length(facG_root_undirected.var)~=n)
        p = rand(1,1);
        G_undirected=random_graph(n,p);
        facG_root_undirected = adjG2facG(G_undirected,min_dim,max_dim);
    end
    while(length(facG_root_directed.var)~=n)
        p = rand(1,1);
        G_directed=random_directed_graph(n,p);
        facG_root_directed = adjG2facG(G_directed,min_dim,max_dim);
    end
    Stats_undirected(1,i)=numedges(G_undirected);
    Stats_undirected(2,i)=num_conn_comp(G_undirected);
    Stats_undirected(3,i)=diameter(G_undirected);
    Stats_undirected(4,i)=ave_path_length(G_undirected);
    Stats_undirected(5,i)=graph_radius(G_undirected);
    Stats_undirected(6,i)=loops3(G_undirected);
    Stats_undirected(7,i)=length(loops4(G_undirected));
    Stats_undirected(8,i)=num_loops(G_undirected);
    Stats_undirected(9,i)=mean(cardi(facG_root_undirected));
    Stats_undirected(10,i)=max(cardi(facG_root_undirected));
    
    Stats_directed(1,i)=numedges(G_directed);
    Stats_directed(2,i)=num_conn_comp(G_directed);
    Stats_directed(3,i)=diameter(G_directed);
    Stats_directed(4,i)=ave_path_length(G_directed);
    Stats_directed(5,i)=graph_radius(G_directed);
    Stats_directed(6,i)=loops3(G_directed);
    Stats_directed(7,i)=length(loops4(G_directed));
    Stats_directed(8,i)=num_loops(G_directed);
    Stats_directed(9,i)=mean(cardi(facG_root_directed));
    Stats_directed(10,i)=max(cardi(facG_root_directed));
    
    
    B_PR_undirected = [];
    I_PR_undirected = [];
    ifconverged_PR_undirected = [];
    B_PR_directed = [];
    I_PR_directed = [];
    ifconverged_PR_directed = [];
    for prior_i=1:length(prior)
        facG_undirected=initialize(facG_root_undirected,prior(prior_i));
        [facG_undirected, iters_undirected,~] = run_loopy_bp_parallel(facG_undirected, max_iters, conv_tol, false);
        [nodeMarg_par_undirected] = get_beliefs(facG_undirected);
        B_PR_undirected = [B_PR_undirected nodeMarg_par_undirected];
        I_PR_undirected = [I_PR_undirected iters_undirected];
        ifconverged_PR_undirected = [ifconverged_PR_undirected iters_undirected<max_iters];
        facG_directed=initialize(facG_root_directed,prior(prior_i));
        [facG_directed, iters_directed,~] = run_loopy_bp_parallel(facG_directed, max_iters, conv_tol, false);
        [nodeMarg_par_directed] = get_beliefs(facG_directed);
        B_PR_directed = [B_PR_directed nodeMarg_par_directed];
        I_PR_directed = [I_PR_directed iters_directed];
        ifconverged_PR_directed = [ifconverged_PR_directed iters_directed<max_iters];
    end
    flag_undirected = 0;
    for j=1:size(B_PR_undirected,2)-1
        if(belief_diff(B_PR_undirected(:,j),B_PR_undirected(:,j+1))>conv_tol)
            flag_undirected = 1;
        end
    end
    if(flag_undirected ==0)
        I_undirected = [I_undirected sum(I_PR_undirected(ifconverged_PR_undirected==1))/length(I_PR_undirected(ifconverged_PR_undirected==1))];
        ifconverged_undirected = [ifconverged_undirected;ifconverged_PR_undirected];
    end
    valid_undirected = [valid_undirected flag_undirected==0];
    
    
    flag_directed = 0;
    for j=1:size(B_PR_directed,2)-1
        if(belief_diff(B_PR_directed(:,j),B_PR_directed(:,j+1))>conv_tol)
            flag_directed = 1;
        end
    end
    if(flag_directed ==0)
        I_directed = [I_directed sum(I_PR_directed(ifconverged_PR_directed==1))/length(I_PR_directed(ifconverged_PR_directed==1))];
        ifconverged_directed = [ifconverged_directed;ifconverged_PR_directed];
    end
    valid_directed = [valid_directed flag_directed==0];
end
    Stats_undirected(Stats_undirected==inf)=0;
    Stats_directed(Stats_directed==inf)=0;
    
for i=1:length(Stats_name)
    stat_undirected = Stats_undirected(i,:);
    stat_directed = Stats_directed(i,:);
    [stat_undirected_x, I_undirected_y]=ave_scatter(stat_undirected(valid_undirected==1),I_undirected(valid_undirected==1));
    [stat_directed_x, I_directed_y]=ave_scatter(stat_directed(valid_directed==1),I_directed(valid_directed==1));
    figure
    scatter(stat_undirected_x,I_undirected_y,15,'b','filled');
    hold on
    scatter(stat_directed_x,I_directed_y,15,'r','filled');
    hold off
    grid on
    legend('undirected','directed')
    xlabel(Stats_name{i})
    ylabel('# of iterations')
    title(sprintf('Convergence rate analysis n=%d',n))
    lsline
end

for i=1:length(Stats_name)
    stat_undirected = Stats_undirected(i,:);
    stat_directed = Stats_directed(i,:);
%% Averaged Version Plot   
    [stat_undirected_x, I_undirected_y]=ave_scatter(stat_undirected(valid_undirected==1),I_undirected(valid_undirected==1));
    [stat_directed_x, I_directed_y]=ave_scatter(stat_directed(valid_directed==1),I_directed(valid_directed==1));
    figure
    scatter(stat_undirected_x,I_undirected_y,15,'b','filled');
    hold on
    scatter(stat_directed_x,I_directed_y,15,'r','filled');
    hold off
    grid on
    legend('undirected','directed')
    xlabel(Stats_name{i})
    ylabel('# of iterations')
    title(sprintf('Convergence rate analysis n=%d',n))
    lsline
%% Original Version Plot   
    figure
    scatter(stat_undirected(valid_undirected==1),I_undirected(valid_undirected==1),15,'b','filled');
    hold on
    scatter(stat_directed(valid_directed==1),I_directed(valid_directed==1),15,'r','filled');
    hold off
    grid on
    legend('undirected','directed')
    xlabel(Stats_name{i})
    ylabel('# of iterations')
    title(sprintf('Convergence rate analysis n=%d',n))
end

save('Stat.mat','Stats_name','Stats_undirected','Stats_directed','I_undirected','I_directed','valid_undirected','valid_directed','n')
converge_rate_wrt_prior_undirected = sum(ifconverged_undirected,1)/size(ifconverged_undirected,1);
converge_rate_wrt_prior_directed = sum(ifconverged_directed,1)/size(ifconverged_directed,1);
