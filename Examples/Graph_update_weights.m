function G = Graph_update_weights(G,zero_threshold)
% for undirected simple graph
% zero_threshold 

n = nnz(G.W)/2;

w=(rand(n,1)*(1-zero_threshold))+zero_threshold;
B=triu(full(G.W));
k=1;

N = size(G.W,1);

for i=1:size(N,1)
    for j=i+1:size(N,2)
        if B(i,j)>0
            B(i,j)=w(k);
            k=k+1;
        end
    end
end


B=B+B';
Wnew = sparse(B);
G.W=Wnew;
G = gsp_graph_default_parameters(G);
end
