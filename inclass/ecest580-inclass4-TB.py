# In Class Assignment 4
# Tyler Bradley
# 2/1/2019

class TreeNode:
    def __init__(self, val):
        self.left = None
        self.right = None
        self.val = val
        
def TreeGenerator(A):
    if A is None or len(A) == 0:
        return None
    mid = len(A) // 2
    root = TreeNode(A[mid])
    root.left = TreeGenerator(A[:mid])
    root.right = TreeGenerator(A[mid+1:])
    return root

def go_left(tree):
    left_tree = tree.left
    left_tree.top = tree
    return left_tree

def go_right(tree):
    right_tree = tree.right
    right_tree.top = tree
    return right_tree


A = [1, 2, 3, 4, 5, 6, 7]

tree = TreeGenerator(A)


output = []
for i in A:
    current_tree = tree
    while len(output) < i:
        if i < current_tree.val:
            current_tree = go_left(current_tree)
        elif i > current_tree.val:
            current_tree = go_right(current_tree)
        else:
            output.append(current_tree.val)
            

output
# Out[57]: [1, 2, 3, 4, 5, 6, 7]     
        
        