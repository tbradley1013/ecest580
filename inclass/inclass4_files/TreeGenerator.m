function root = TreeGenerator(input_list)
%% create a BST from a sorted list
    if isempty(input_list)
        root = [];
        return
    end
    mid = ceil(length(input_list) / 2);
    root = TreeNode;
    root.Val = input_list(mid);
    root.Left = TreeGenerator(input_list(1:mid - 1));
    root.Right = TreeGenerator(input_list(mid + 1:end));
    return
end
