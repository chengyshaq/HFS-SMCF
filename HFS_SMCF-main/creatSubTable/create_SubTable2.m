%% creatSubTable创建树模型
% Written by Yu Wang
% Modified by Hong Zhao
% Modified by Haoyang Liu
%% Creat subtable
function [DataMod,LabelMod,IndxMod]=create_SubTablezh2(dataset, tree)
Data = dataset(:,1:end);%样本特征矩阵
Label =  dataset(:,end);%样本标签矩阵
[numTrain,~] = size(dataset);%样本数
internalNodes = tree_InternalNodes(tree);%找树的中间节点
internalNodes(find(internalNodes==-1))=[];%中间节点为-1的置空
indexRoot = tree_Root(tree);% The root of the tree树的根节点
noLeafNode =[internalNodes;indexRoot];%非叶子节点
for i = 1:length(noLeafNode)%循环遍历每一个非叶子节点
    cur_descendants = tree_Descendant(tree, noLeafNode(i));%第i个非叶子节点的子节点
    ind_d = 1;  % index for id subscript increment  id下标增量的索引
    id = [];        % data whose labels belong to the descendants of the current nodes其标签属于当前节点的后代的数据
    for n = 1:numTrain  %循环遍历每一个样本
        if (ismember(Label(n), cur_descendants) ~= 0)%判断该样本是否属于第i个节点
            id(ind_d) =  n;
            ind_d = ind_d +1;
        end
    end
    Label_Uni_Sel = Label(id,:);%选出第i个中间节点的样本对应的标签
    DataSel = Data(id,:);     %select relative training data for the current classifier为当前分类器选择相应的训练数据

    numTrainSel = size(Label_Uni_Sel,1);%第i个非叶子节点对应的样本数
    LabelUniSelMod = mylabel_modify_MLNP(Label_Uni_Sel, noLeafNode(i), tree);%% 得到标签树模型
    % Get the sub-training set containing only relative nodes得到只包含相对节点的子训练集
    ind_tdm = 1;
    index = [];     % data whose labels belong to the children of the current nodes标签属于当前节点的子节点的数据
    children_set = get_children_set(tree, noLeafNode(i));%得到当前节点的叶子节点
    for ns = 1:numTrainSel%循环遍历每一个实例
        if (ismember(LabelUniSelMod(ns), children_set) ~= 0)
            index(ind_tdm) =  ns;
            ind_tdm = ind_tdm +1;%统计实例个数
        end
    end
    DataMod{noLeafNode(i)} = DataSel(index, :);   % Find the sub training set of relative to-be-classified nodes
                                                    %找到相对的分类节点的子训练集
    LabelMod{noLeafNode(i)} = LabelUniSelMod(index, :);%叶子节点标签模型
    IndxMod{noLeafNode(i)} = id;%属于第i个非叶子节点的样本索引模型
end
end