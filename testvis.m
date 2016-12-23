clear;
max_iters = 500; conv_tol = 1e-6; total_iters = 500;
min_dim = 2; max_dim = 5;
range_dim = max_dim-min_dim;
n = 10;
nodeMarg_undirected = {};
nodeMarg_directed = {};
G_undirected_obj = {};
G_directed_obj ={};
for i=1:length(n)
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
    G_undirected_obj = [G_undirected_obj G_undirected];
    G_directed_obj = [G_directed_obj G_directed];
    facG_undirected=initialize(facG_root_undirected,1);
    facG_directed=initialize(facG_root_directed,1);
    [facG_undirected, iters_undirected, nodeMarg_par_undirected] = run_loopy_bp_parallel(facG_undirected, max_iters, conv_tol, true);
    [facG_directed, iters_directed, nodeMarg_par_directed] = run_loopy_bp_parallel(facG_directed, max_iters, conv_tol, true);
    if(~(iters_undirected<max_iters&&iters_directed<max_iters))
        disp('fail to converge, retry pls');
        break;
    end
    C_undirected=[];
    C_directed=[];
    dim_undirected=zeros(1,size(nodeMarg_par_undirected,1));
    dim_directed=zeros(1,size(nodeMarg_par_directed,1));
    for j=1:size(nodeMarg_par_undirected,1)
        color_undirected = [];
        for m=1:size(nodeMarg_par_undirected,2)
            tmp_undirected=[1:length(nodeMarg_par_undirected{j,m})] * nodeMarg_par_undirected{j,m};
            color_undirected = [color_undirected tmp_undirected];
        end
        dim_undirected(1,j) = length(nodeMarg_par_undirected{j,1});
        C_undirected=[C_undirected;color_undirected];
    end
    for q=1:size(nodeMarg_par_directed,1)
        color_directed = [];
        for k=1:size(nodeMarg_par_directed,2)
            tmp_directed=[1:length(nodeMarg_par_directed{q,k})] * nodeMarg_par_directed{q,k};
            color_directed = [color_directed tmp_directed];
        end
        dim_directed(1,q) = length(nodeMarg_par_directed{q,1});
        C_directed=[C_directed;color_directed];
    end
    marg_p_undirected = [];
    marg_p_directed = [];
    [~,max_dim_idx_undirected] = max(dim_undirected);
    [~,max_dim_idx_directed] = max(dim_directed);
    for u=1:size(nodeMarg_par_undirected,2)
        marg_p_undirected = [marg_p_undirected nodeMarg_par_undirected{max_dim_idx_undirected,u}];
    end
    for y=1:size(nodeMarg_par_directed,2)
        marg_p_directed = [marg_p_directed nodeMarg_par_directed{max_dim_idx_directed,y}];
    end
    figure
    subplot(2,1,1)
    plot(marg_p_undirected')
    xlabel('# of iterations')
    ylabel('Marginal Probability')
    title('Convergence Visulization(undirected)')
    subplot(2,1,2)
    plot(marg_p_directed')
    xlabel('# of iterations')
    ylabel('Marginal Probability')
    title('Convergence Visulization(directed)')
    
    figure
    subplot(2,1,1)
    plot(C_undirected')
    xlabel('# of iterations')
    ylabel('Marginal Expectation')
    title('Convergence Visulization(undirected)')
    subplot(2,1,2)
    plot(C_directed')
    xlabel('# of iterations')
    ylabel('Marginal Expectation')
    title('Convergence Visulization(directed)') 
end


