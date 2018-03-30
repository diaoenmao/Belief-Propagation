function [G, iters, node_marg] = run_loopy_bp_parallel( G, max_iters, conv_tol,ifvis )
  node_marg =[];
  for iters = 1:max_iters
    for i = 1:length(G.fac)
      G.fac(i).oldoutgoing=G.fac(i).outgoing;
      mnum = length(G.fac(i).incoming);
      for j = 1:mnum
          nextShape = size(G.fac(i).p);
          nextShape(j)=[];
          nextShape=[1 nextShape];
          tmp=repmat(G.fac(i).incoming{j},nextShape);
          G.fac(i).incoming{j}=ipermute(tmp,[j [1:j-1,j+1:length(nextShape)]]);              
      end
      for j = 1:mnum
          curr = G.fac(i).incoming;
          curr(j)=[];
          newmessage=G.fac(i).p;
          for m = 1:length(curr)
            newmessage = newmessage .* curr{m};
          end
          tmp=permute(newmessage,[j,[1:j-1,j+1:length(size(newmessage))]]);
          dim = 2:length(size(newmessage));
          for n = 1:length(dim)
            tmp=sum(tmp,dim(n));
          end
          G.fac(i).outgoing{j}=tmp/sum(tmp);
      end     
    end
    for i=1:length(G.fac)
        for j=1:length(G.fac(i).outgoing)
            varidx=G.fac(i).nbrs_var(j);
            G.var(varidx).incoming{G.var(varidx).nbrs_fac==i}=G.fac(i).outgoing{j};
        end
    end
    for i=1:length(G.var)
        if(length(G.var(i).nbrs_fac)>1)
            G.var(i).oldoutgoing=G.var(i).outgoing;
            for j=1:length(G.var(i).incoming)
                curr = G.var(i).incoming;
                curr(j)=[];
                newmessage=1;
                for m = 1:length(curr)
                    newmessage = newmessage .* curr{m};
                end
                G.var(i).outgoing{j}=newmessage/sum(newmessage);
            end
        end
    end
    for i=1:length(G.var)
        for j=1:length(G.var(i).outgoing)
            facidx=G.var(i).nbrs_fac(j);
            G.fac(facidx).incoming{G.fac(facidx).nbrs_var==i}=G.var(i).outgoing{j};
        end
    end
    delta=[];
    for i=1:length(G.var)
        for j=1:length(G.var(i).oldoutgoing)
            delta=[delta sum(abs(G.var(i).oldoutgoing{j}-G.var(i).outgoing{j}))];
        end
    end
    for i=1:length(G.fac)
        for j=1:length(G.fac(i).oldoutgoing)
            delta=[delta sum(abs(G.fac(i).oldoutgoing{j}-G.fac(i).outgoing{j}))];
        end
    end
    if(ifvis)
        node_marg = [node_marg get_beliefs(G)];
    end
    if(all(delta<conv_tol))
        break;
    end  
  end
end

