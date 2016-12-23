function out = cardi(facG)
facN = length(facG.fac);
out = zeros(1,facN);
for i=1:facN
    out(i) = length(size(facG.fac(i).p));   
end
end