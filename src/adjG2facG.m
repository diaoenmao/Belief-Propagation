function facG = adjG2facG(G,mindim,maxdim)
if(issymmetric(G))
    %Undirected Graph to MRF(pairwise factor graph)
    G = triu(G,1);
    facG = init_graph();
    varN = [];
    for i=1:size(G,1)
        for k=i:size(G,2)
            if(G(i,k)==1)
                if(mindim~=maxdim)
                    dims = randi(maxdim-(mindim-1),1,2)+(mindim-1);
                else
                    dims = [mindim mindim];
                end
                ids=zeros(1,2);
                if(i~=k)
                    if(any(varN==i))
                        varid = varN==i;
                        ids(1) = facG.var(varid).id;
                        dims(1) = facG.var(varid).dim;
                    else
                        varN = [varN i];
                        [facG,ids(1)]=add_varnode(facG,strcat('v',int2str(i)),dims(1));
                    end
                    if(any(varN==k))
                        varid = varN==k;
                        ids(2) = facG.var(varid).id;
                        dims(2) = facG.var(varid).dim;
                    else
                        varN = [varN k];
                        [facG,ids(2)]=add_varnode(facG,strcat('v',int2str(k)),dims(2));
                    end
                    p=rand(dims(2),dims(1));
                    facG=add_facnode(facG,p,ids(2),ids(1));
                else
                    if(any(varN==i))
                        varid = varN==i;
                        ids(1) = facG.var(varid).id;
                        dims(1) = facG.var(varid).dim;
                    else
                        varN = [varN i];
                        [facG,ids(1)]=add_varnode(facG,strcat('v',int2str(i)),dims(1));
                    end
                        p=rand(dims(1),dims(1));
                        facG=add_facnode(facG,p,ids(1),ids(1));
                end
            end
        end
    end
else
    %Directed Graph to hiorder potential(general factor graph)
    facG = init_graph();
    varN = [];
    pointed = [];
    for i=1:size(G,1)
        for k=1:size(G,2)
           if(G(i,k)==1)
               if(mindim~=maxdim)
                   dims = randi(maxdim-(mindim-1),1,2)+(mindim-1);
               else
                   dims = [mindim mindim];
               end
               ids=zeros(1,2);
               if(any(pointed==k))
                   facid = find(pointed==k);
                   if(any(varN==i))
                       varid = varN==i;
                       ids(1) = facG.var(varid).id;
                       facG.var(ids(1)).nbrs_fac = [facG.var(ids(1)).nbrs_fac;facid];
                   else
                      varN = [varN i];
                      if(i~=k)
                        [facG,ids(1)]=add_varnode(facG,strcat('v',int2str(i)),dims(1));
                      end
                      facG.var(ids(1)).nbrs_fac = facid;
                   end
                   facG.fac(facid).nbrs_var = [facG.fac(facid).nbrs_var ids(1)];
                   varids = facG.fac(facid).nbrs_var;
                   dims = [];
                   for j=1:length(varids)
                       dims=[dims facG.var(varids(j)).dim];
                   end
                   facG.fac(facid).p = rand(dims);
               else
                   pointed=[pointed k];
                   if(i~=k)
                       if(any(varN==i))
                           varid = varN==i;
                           ids(1) = facG.var(varid).id;
                           dims(1) = facG.var(varid).dim;
                       else
                           varN = [varN i];
                           [facG,ids(1)]=add_varnode(facG,strcat('v',int2str(i)),dims(1));
                       end
                       if(any(varN==k))
                           varid = varN==k;
                           ids(2) = facG.var(varid).id;
                           dims(2) = facG.var(varid).dim;
                       else
                           varN = [varN k];
                           [facG,ids(2)]=add_varnode(facG,strcat('v',int2str(k)),dims(2));
                       end
                       p=rand(dims(2),dims(1));
                       facG=add_facnode(facG,p,ids(2),ids(1));
                   else
                       if(any(varN==i))
                           varid = varN==i;
                           ids(1) = facG.var(varid).id;
                           dims(1) = facG.var(varid).dim;
                       else
                           varN = [varN i];
                           [facG,ids(1)]=add_varnode(facG,strcat('v',int2str(i)),dims(1));
                       end
                       p=rand(dims(1),dims(1));
                       facG=add_facnode(facG,p,ids(1),ids(1));
                   end
               end
           end            
        end
    end
    non_pointed = varN(~ismember(varN,pointed));
    for i=1:length(non_pointed)
        varid = varN==non_pointed(i);
        p=rand(facG.var(varid).dim,1);
        facG = add_facnode(facG,p,facG.var(varid).id);
    end
end
end