#ifndef _TREE_H
#define _TREE_H

#include <iostream>
#include <vector>
using namespace std;

namespace util {

enum NodeType {
  inner,
  leaf
};

class TreeNode {
  public:
  TreeNode* father;
  TreeNode* child;
  TreeNode *prevSib, *nextSib;
  NodeType type;
  int data;   
  TreeNode() {father = NULL; child = NULL; prevSib = NULL; nextSib = NULL; }
  TreeNode(int i, NodeType t) {
    data = i; 
    type = t;
    father = NULL;
    child = NULL;
    prevSib = NULL;
    nextSib = NULL;
  }
  TreeNode(int i) {
    data = i; 
    type = leaf;
    father = NULL;
    child = NULL;
    prevSib = NULL;
    nextSib = NULL;
  }
  TreeNode* addChild(int i) {
    TreeNode* tn = new TreeNode(i); 
    // loop over child, add the new node into its children list
    TreeNode* ch = child;
    // if no children yet, the added one is its first child
    if (child==NULL) child = tn; 
    else { 
      while (ch->nextSib!=NULL) ch = ch->nextSib;
      ch->nextSib = tn;
      tn->prevSib = ch;
    }
    // set father pointer
    tn->father = this; 
    this->type = inner;
    return tn;
  }
};

class Tree {
  public:
  TreeNode* root;
  Tree() {}
  Tree(int i) {root = new TreeNode(i);}
  Tree(TreeNode* r) {root = r;}

  void setRoot(TreeNode* r) {root = r;}

  void dfs(TreeNode* r, vector<int>& out) {
    if (r!=NULL) out.push_back(r->data);
    else return;
    TreeNode* ch = r->child;
    while (ch!=NULL) {
      dfs(ch,out);
      ch = ch->nextSib;
    }
  }

  void find(TreeNode* r, int i, TreeNode*& f) {
    if (r->data==i) f = r;
    else {
      TreeNode* ch = r->child;
      while (ch!=NULL) {
        find(ch,i,f);
        ch = ch->nextSib;
      }
    }
  } 
  
  void del(TreeNode* t) {
    if (t==NULL) return;
    if (t->father==NULL) return; // not allowed
    if (t->prevSib==NULL) {
      t->father->child = t->nextSib;
      if (t->nextSib!=NULL) t->nextSib->prevSib = NULL;
    }
    else {
      t->prevSib->nextSib = t->nextSib;
    }
  }
};
}

#endif
