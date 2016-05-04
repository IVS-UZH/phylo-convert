# Function for inserting a new tip into a phylo tree by splitting off the hierarchy at a specific age (cumulative edge length)
#
# Written by Taras Zakharko based on specs by Balthasar Bickel, [2016-04-26]
# 
# Specifically, given a hierarchy 
#
#  N0 -- N1 -- N2 -- N3 -- reference.tip  
#
# The function will insert (if nessesary) a new split point (node) into the hierarchy and branch the new.tip off it
#
#
# N0 -- N1 -- N' -- N2 -- N3 -- reference.tip 
#              \
#               \_____ new.tip
#
# such that length(N', reference.tip) = split.branch.length and length(N', new.tip) = new.tip.branch.length
# The length of the branches are preserved, i.e. original length(N1, N2) == length(N1, N') + length(N', N2)
 
graft <- function(
  tree,                 # a phylo object
  new.tip,              # the name of the new tip to be inserted
  reference.tip,        # an existing tip in whose hierarchy the new tip should be inserted
  split.branch.length,  # the cumulative edge length (age) from the reference.tip to the new.tip ancestor, i.e. how many years before the reference.tip did the new tip branch off?
  new.tip.branch.length # the length of the edge between new.tip and its parent, i.e. how long did the new tip evolved on its own?
  ) {
    # validate the input
    stopifnot(inherits(tree, "phylo"))
    stopifnot(!is.null(tree$edge.length))
    stopifnot(reference.tip %in% tree$tip.label)
    stopifnot(inherits(new.tip, 'character'))
    stopifnot(!new.tip %in% tree$tip.label)
    stopifnot(inherits(split.branch.length, 'numeric'))
    stopifnot(inherits(new.tip.branch.length, 'numeric'))
    
    # a utility function for retrieving the parent of a node or NA if no such parent exist
    getParent <- function(node) {
      id <- tree$edge[tree$edge[, 2] %in% node, 1]
      if(length(id)!=1) NA else id
    }
    # a utility function for retrieving edge id. Precondition: edge exists
    getEdgeID <- function(parent, child) { 
      id <- which((tree$edge[, 1] %in% parent) & (tree$edge[, 2] %in% child))
      stopifnot(length(id) == 1)
      id
    }
    
    # a utility function for computing the total branch length from a node to tip
    # precondition: node is an ancestor of tip
    blength <- function(from, to) {
      # walk the branches up until we hit the node we seek
      node <- to
      cumulative.length <- 0
      
      while(node != from) {
        parent <- getParent(node)
        if(is.na(parent)) stop("Unexpectedly reached root when traversing the tree")
          
        # increase the length by the branch length   
        cumulative.length <-  cumulative.length + tree$edge.length[getEdgeID(parent, node)]
        
        # and walk one step up
        node <- parent 
      }
      
      return(cumulative.length)
    }
    
    
    
    # ------------------- algorithm start
    
    # 1. insert the new tip and adjust the node ids accordingly (so that we don't mess with them later)
    #    Specifically, the node ids need to be moved to the right (+1)  
    tree$edge <- with(tree, ifelse(edge > length(tip.label), edge+1, edge))
    tree$tip.label <- c(tree$tip.label, new.tip)
    
    # get the tip ids 
    reference.tip.id <- which(tree$tip.label %in% reference.tip)
    new.tip.id      <- which(tree$tip.label %in% new.tip)
    
    # 2. locate the nodes A and B between which the new split should be inserted
    #    this is the node such that blength(B, reference.tip) < split.branch.length and 
    #    blength(A, reference.tip) >= split.branch.length
    # 
    #    edge cases: the split point coincides with A (do not insert a new node in this case)
    #                B is the tree root (insert a new node as a new root and branch new tip off it)
    
    # init B to the reference.tip 
    B <- reference.tip.id
    # and a to its ancestor
    A <- getParent(B)
        
    # we go up until our search conditions are satisfied
    while(!(is.na(A) || # B has no ancestor 
            (blength(B, reference.tip.id) < split.branch.length && # split point is between A and B
             blength(A, reference.tip.id) >= split.branch.length)
          )) {
              # go one step up
              B <- A
              A <- getParent(B)
           }
    
    # edge case 1: B is root (A is NA)
    if(is.na(A)) {
      # we insert a new node N as the new root and branch both B and new tip off it
      # the length of the A -- B branch is split.branch.length - blength(B, reference.tip.id)
      tree$Nnode <- tree$Nnode + 1
      N <- length(tree$tip.label) + tree$Nnode
      
      # branch N -- B
      tree$edge        <- rbind(tree$edge, c(N, B))
      tree$edge.length <- c(tree$edge.length, split.branch.length - blength(B, reference.tip.id))
      
      # branch N -- new.tip
      tree$edge        <- rbind(tree$edge, c(N, new.tip.id))
      tree$edge.length <- c(tree$edge.length, new.tip.branch.length)
    } else
    # edge case 2: split point is exactly A
    if(blength(A, reference.tip.id) == split.branch.length) {
      # we simply branch new.tip directly off A
      tree$edge        <- rbind(tree$edge, c(A, new.tip.id))
      tree$edge.length <- c(tree$edge.length, new.tip.branch.length)
    } else {
      # we need to insert a new node between A and B
      tree$Nnode <- tree$Nnode + 1
      N <- length(tree$tip.label) + tree$Nnode
      
      # compute the new branch lengths
      A.N.length <- blength(A, reference.tip.id) - split.branch.length
      N.B.length <- split.branch.length - blength(B, reference.tip.id)
      
      # delete the old branch
      tree$edge.length <- tree$edge.length[-getEdgeID(A, B)]
      tree$edge        <- tree$edge[-getEdgeID(A, B), ]
      
      # branch A -- N 
      tree$edge        <- rbind(tree$edge, c(A, N))
      tree$edge.length <- c(tree$edge.length, A.N.length)
      
      # branch N -- B
      tree$edge        <- rbind(tree$edge, c(N, B))
      tree$edge.length <- c(tree$edge.length, N.B.length)
      
      # branch N -- new.tip
      tree$edge        <- rbind(tree$edge, c(N, new.tip.id))
      tree$edge.length <- c(tree$edge.length, new.tip.branch.length)      
    }
    
  
    # return the result
    attr(tree, 'order') <- NULL
    tree
}

# ----- test code here
#
# library(ape)
# library(dplyr)
#
# tt <- read.nexus("mcc.trees")
# tt1 <- read.nexus("ie.tree.full.nex")
#
# extractEdges <- function(phylo) {
#   edges <- data.frame(as.data.frame(phylo$edge), phylo$edge.length)
#   names(edges) <- c('parent', 'child', 'length')
#
#   edges
# }
#
# subset(extractEdges(tt), child %in% c(22, 126))
#
#
# xx <- insert_tip_into_tree_with_ages(tt, 'Middle_English', 'English', 691, 0.9)
# xx <- insert_tip_into_tree_with_ages(xx, 'Middle_Welsh', 'Welsh', 761, 0.1)
# xx <- insert_tip_into_tree_with_ages(xx, 'Middle_Persian', 'Persian', 1577, 0.9)
# xx <- insert_tip_into_tree_with_ages(xx, 'Luwian', 'Hittite', 1000, 1500)
#
#
# subset(extractEdges(xx), child %in% c(166, which(xx$tip.label %in% c('English', 'Middle_English'))))
# subset(extractEdges(xx), child %in% c(167, which(xx$tip.label %in% c('Middle_Welsh', 'Welsh'))))
# subset(extractEdges(xx), child %in% c(102, 99, 169, which(xx$tip.label %in% c('Persian'))))
# subset(extractEdges(xx), child %in% c(169, which(xx$tip.label %in% c('Middle_Persian'))))
# subset(extractEdges(xx), child %in% c(171, which(xx$tip.label %in% c('Luwian', 'Hittite'))))
#
#
#
#
#
#
#

