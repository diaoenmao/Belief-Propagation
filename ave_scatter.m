function [x,y]= ave_scatter(org_x,org_y)
u_x = unique(org_x);
if(length(u_x)==1)
    x=u_x;
    y=mean(org_y(~isnan(org_y)));
    return;
end
dup_idx = u_x(hist(org_x,u_x)>1);
x=[];
y=[];
removed_idx=[];
for i=1:length(org_x)
    if(any(removed_idx==i))
        continue;
    end
    if(any(dup_idx==org_x(i)))
        all_dup_idx = find(org_x==org_x(i));
        other_dup_idx = all_dup_idx(all_dup_idx~=i);
        removed_idx = [removed_idx other_dup_idx];
        tmp = org_y(all_dup_idx);
        ave_dup_y = mean(tmp(~isnan(tmp)));
        x=[x org_x(i)];
        y=[y ave_dup_y];
    else
        x=[x org_x(i)];
        y=[y org_y(i)];
    end
end
end