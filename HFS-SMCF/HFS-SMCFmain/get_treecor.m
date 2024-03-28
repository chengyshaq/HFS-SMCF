function [treecor] = get_treecor(tree,clus,data_array)
%根据树结构计算节点间的相似度,标签相关性
[r,~]=size(tree);
for i=1:r
    n=i;
    while(n~=0)
        label(i,n)=1;
        n=tree(n,1);
    end  
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  net = aabpTrain(data_array(:,1:end-1),lable);
%  preDistribution = aabpPredict(net, data_array(:,1:end-1));
%  label
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:r
    for j=1:r
        treecor(i,j)=corr(label(i,:)',label(j,:)','type','pearson');
    end
end
     a=corr(clus',clus','type','pearson');
     [row, col] = find(isnan(a));  % 找出 A 中所有的 NaN 元素
    for i = 1:length(row)
        a(row(i), col(i)) = 0;  % 将 NaN 元素赋值为 0
    end
     treecor =0.9*treecor+0.1*a;
end

