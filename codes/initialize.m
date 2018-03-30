function G = initialize(G,prior)
for i=1:length(G.fac)
    G.fac(i).outgoing = cell(1,length(G.fac(i).nbrs_var));
    G.fac(i).incoming = cell(1,length(G.fac(i).nbrs_var));
    G.fac(i).oldoutgoing = cell(1,length(G.fac(i).nbrs_var));
    for j = 1:length(G.fac(i).nbrs_var)
        G.fac(i).outgoing{j} = prior*ones(G.var(G.fac(i).nbrs_var(j)).dim,1);
        G.fac(i).incoming{j} = prior*ones(G.var(G.fac(i).nbrs_var(j)).dim,1);
        G.fac(i).oldoutgoing{j} = prior*ones(G.var(G.fac(i).nbrs_var(j)).dim,1);
    end
end
for i=1:length(G.var)
    G.var(i).outgoing = cell(1,length(G.var(i).nbrs_fac));
    G.var(i).incoming = cell(1,length(G.var(i).nbrs_fac));
    G.var(i).oldoutgoing = cell(1,length(G.var(i).nbrs_fac));
    for j = 1:length(G.var(i).nbrs_fac)
        G.var(i).outgoing{j} = prior*ones(G.var(i).dim,1);
        G.var(i).incoming{j} = prior*ones(G.var(i).dim,1);
        G.var(i).oldoutgoing{j} = prior*ones(G.var(i).dim,1);
    end
end
end