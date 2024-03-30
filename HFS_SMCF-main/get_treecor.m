function [treecor] = get_treecor(tree,clus,data_array)
%�������ṹ����ڵ������ƶ�,��ǩ�����
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
     [row, col] = find(isnan(a));  % �ҳ� A �����е� NaN Ԫ��
    for i = 1:length(row)
        a(row(i), col(i)) = 0;  % �� NaN Ԫ�ظ�ֵΪ 0
    end
     treecor =0.9*treecor+0.1*a;
end

