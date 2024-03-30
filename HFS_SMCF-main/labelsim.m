function [labelsimi] = labelsim(tree,data_array)
%根据树结构计算节点间的相似度
[r,~]=size(tree);
label = zeros(size(data_array,1),r);
label_0 = data_array(:,end);
for i=1:size(data_array,1)
    a=label_0(i,1);
    label(i,a)=1;
    b=tree(a,1);
    label(i,b)=1;
    label(i,end)=1;
end
% net = aabpTrain(X2{noLeafNode(i)}(:,1:end-1),Y{noLeafNode(i)});
% preDistribution{noLeafNode(i)} = aabpPredict(net, X2{noLeafNode(i)}(:,1:end-1));

labelsimi=corr(label,label,'type','pearson');
labelsimi(isnan(labelsimi)) = 0;
%%%%%%%%%%%%%%%%%%%%%%
% dis = pdist2(label', label');
% labelsimi = sign(labelsimi);
% labelsimi=labelsimi.*dis;
% labelsimi = mapminmax(labelsimi);
%%%%%%%%%%%%%%%%%%%%%%%

end