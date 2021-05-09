#ifndef AVL_TREE
#define AVL_TREE

#include <inttypes.h>

#define LL 1
#define LR 2
#define RR 3
#define RL 4


struct tree_node {
    uint64_t value;
    unsigned long index;
    char bf;
    char branch; //relative to its parent.
    struct tree_node * parent, * left, * right;
};

struct tree_node * make_avl_tree(struct tree_node ** const, uint64_t, unsigned long *);

#endif
