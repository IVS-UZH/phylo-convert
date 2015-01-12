require(ape)

# this computes the height of the phylo tree
levels.phylo <- function(x)
{
	# initial, very slow but argubly clean version
	# TODO: optimize this
	edges <- x$edge
	roots <- setdiff(1:max(x$edge), edges[, 2])
	
	node.height <- function(node) 
	{
		children <- edges[edges[, 1] == node, 2]
		if(length(children)==0) 
			return(1)
		else
			return(1 + max(sapply(children, node.height)))
	}
	
	max(sapply(roots, node.height))
}


levels.multiPhylo <- function(x) max(sapply(x, levels.phylo))

# x is a phylo object 
# we return a data frame with lowest-level taxon in every row
as.data.frame.phylo <- function(x, row.names = NULL, optional = FALSE, height = levels(x), ...)
{
	# initial, very slow but argubly clean version
	# TODO: optimize this

	# this is our output
	rows <- NULL
	
	# names of the data frame
	df_names <- paste0('Level', rev(1:height))
	
	# the nodes (I don't differentiate between parents and leafs here)
	nodes <- c(x$tip.label, x$node.label)
	
	# the edge matrix, we use this to walk the tree
	edges <- x$edge
	
	# this function recursively walks the tree until its at a leaf
	tree_walker <- function(path)
	{
		# cat(paste0('walking...', deparse(path), '\n'))
		
		
		# get the child nodes
		last_node = path[length(path)]
		children <- edges[edges[, 1] == last_node, 2]
		
		if (length(children) == 0)
		{	
			# cat('we are here', deparse(path), '\n')
			
			# we have a full path!
			# now we want to expand it using NAs if its not the full length
			if(length(path)<height) 
				p1 <- nodes[c(path[-length(path)], rep(NA, height - length(path)), last_node)]
			else
				p1 <- nodes[path]
			
			# print(p1)
			
			p1 <- as.data.frame(as.list(p1))
			names(p1) <- df_names
			# cat('Output line:\n')
			# print(p1)
			
			rows <<- rbind(rows, p1)
		}
		else for(child in children)
		{
			p1 <- c(path, child)
			tree_walker(p1)
		}
	}
		
	# find the root nodes, these are the nodes which don't have any parents
	roots <- setdiff(seq_along(nodes), edges[, 2])

	# and walk them
	for(root in roots) tree_walker(root)
		
	rownames(rows) <- NULL
		
	return(rows)	
}

as.data.frame.multiPhylo <- function(x, row.names = NULL, optional = FALSE, height = levels(x), ...)
{
	# # test if multicore is installed and available
	# if(is.element('parallel', installed.packages()[,1])) 
	# {
	# 	require(parallel)
	# 	applyfunc <- function(x, f) mclapply(x, f)
	# }
	# else applyfunc <- function(x, f) lapply(x, f)
	# 
	# do.call(rbind, applyfunc(x, function(tree) as.data.frame.phylo(tree, height = height)))
	do.call(rbind, lapply(x, as.data.frame.phylo, height = height))
	
}



propagate.node.to.tip <- function (x, node = NULL, filter = NULL, ...)  UseMethod("propagate.node.to.tip")


propagate.node.to.tip.phylo <- function(x, root = NULL, filter = NULL, ...)
{
	# some checking + identifying the node
	if(is.null(filter) && is.null(root))
	{
		stop('drop.root requires either a list of roots or a filter function to proceed!')
	} 
	else if(!is.null(filter))
		root.id <- which(sapply(x$node.label, filter))
	else
		root.id <- which(x$node.label %in% root)
	
	# if no hit, just quit
	if(length(root.id)==0) return(x)
		
	# some counts
	old_tip_count <- length(x$tip.label)
	new_tip_count <- length(x$tip.label) + length(root.id)
	node_offset <- length(root.id)
	
	# what we need to do
	# a) offset all node ids by count(new_tips) - count(old_tips)
	# the node ids are identified as being under the old tip id
	x$edge <- ifelse(x$edge > old_tip_count, x$edge + node_offset, x$edge)	
	# b) add the new tips
	x$tip.label <- c(x$tip.label, x$node.label[root.id])	
	# c) add the new edges
	for(i in seq_along(root.id))
		x$edge <- rbind(x$edge, c(root.id[i] + new_tip_count, i + old_tip_count))
	# and re-sort it!
	x$edge <- x$edge[order(x$edge[, 1], x$edge[, 2]), ]
	
	# d) add the edge length
	x$edge.length <- c(x$edge.length, rep(1, length(root.id)))
		
	attr(x, 'order') <- NULL
	
	x
}

propagate.node.to.tip.multiPhylo <- function(x, node = NULL, filter = NULL, ...)
{
	structure(lapply(x, propagate.node.to.tip.phylo, node = node, filter = filter), class = 'multiPhylo')
}

# Converts unrolled genealogical representations to a tree representation (a phylo object)
#
#   
#
as.phylo.data.frame <- function(x, levels = names(x), ...)
{
	# take only the information about the hierarchies
	x <- unique(x[, levels, drop = F])
	
	# make sure there are no factors
	for(l in levels)
		x[[l]] <- as.character(l)
	
	
	# it is possible that some taxa names are duplicated within a hierarchy
	# we need to disambiguate those
	ambiguous_names <- na.omit(unique(unlist(apply(x, 1, function(row)  row[duplicated(row)]))))
	for(n in ambiguous_names)
	{
		# where does that name occur?
		n_levels <-  levels[apply(x == n, 2, any, na.rm=T)]
		warning(paste0("Ambiguous label '", n, "' replaced by ", paste0("'", n, " ", n_levels, "'", collapse=',')))
		for(l in n_levels) 
			x[[l]] <- ifelse(x[[l]] %in% n, paste(n, l), as.character(x[[l]]))
	}

			
	# split the data frame into trees (with unique roots) and convert those to phylo
	trees <- lapply(split(x, x[, 1], drop = T), .data.frame.to.phylo)
	
	# we either return a phylo or a multiPhylo depending on how many roots we have
	if(length(trees) == 1)
		trees[[1]]
	else
		structure(trees, class = 'multiPhylo')
		 
	
}

.data.frame.to.phylo <- function(x)
{
	# first step is to convert the unrolled (flat) hierarchies into direct dominance pairs (graph edges)
	edges <- NULL
	
		
	for(row in 1:nrow(x))
	{
		# init the current parrent to NA
		parent <- NA
		
		# run thourhg the hierarchy
		for(col in 1:ncol(x))
		{
			item <- as.character(x[row, col])
			
			# only continue if the current item is not an NA
			if(!is.na(item))
			{
				# add a new edge if required
				if(!is.na(parent))
					edges <- rbind(edges, c(parent, item))
				
				# the current item becomes the new parent
				parent <- item	
			}
		}
	}

	
	edges <- unique(edges)
	edges <- edges[!(edges[, 1] == edges[, 2]), ]
		
    
	# split the nodes into tips and non-terminal nodes for phylo format
	# tips are nodes which are never parents
	tips <- unique(setdiff(edges[, 2], edges[, 1]))
	
	# nodes are nodes which are not tips
	nodes <- unique(setdiff(edges, tips))
	
	# the id of tips come before the ids of nodes
	node_names <- c(tips, nodes)

	# and now recode the edges to ids
	edges <- matrix(match(edges, node_names), ncol=2)
	

	structure(list(
		edge = edges, 
		Nnode = length(nodes), 
		tip.label = tips, 
		edge.length = rep(1.0, nrow(edges)), 
		node.label = nodes), 
	class = 'phylo')
}

