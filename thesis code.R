################################################################################
# Libraries
library(Matrix)
library(rhdf5)
library(tidyverse)
library(igraph)
library(tidygraph)
library(ggthemes)
library(ggraph)
library(ergm)
library(intergraph)
#install.packages("NetSwan")
#---------------------------

out.h5ls <- h5ls("cons_locs_pathways_mc4_Column.h5")
str(out.h5ls)
head(out.h5ls)

#estrazione informazioni relative al connettoma
id <- out.h5ls$dim != "" 
cMat.info <- out.h5ls[id, ]

# individuazione del gruppo "from" e gruppi "to
strsplit(cMat.info$group, split = "/") |>
  lapply(\(x) x[3:4]) -> group.info
do.call(rbind, group.info) -> group.info
colnames(group.info) <- c("group_to", "group_from")
#cMat.info <- cMat.info[, -1L]
cMat.info <- cbind(group.info[, 2:1], cMat.info)
cMat.info <- as_tibble(cMat.info)

cMat.info

cMat.info |>
  pull(group_to) |>
  unique() -> groupTo


# extraction of information relating to the sub-population 1
i <-  1
cMat.info |>
  filter(group_to == groupTo[1] & is.na(group_from)) |>
  print(n = Inf)

 

cMat.info |>
  filter(group_to == groupTo[3] & !is.na(group_from)) |>
  print(n = Inf)

Adj <- vector(mode = "list", length = length(groupTo))
names(Adj) <- groupTo
for (i in seq_along(groupTo)) {
  cMat.info |>
    filter(group_to == groupTo[i] & !is.na(group_from)) |>
    pull(group) -> group_path
  mapply(\(name) h5read(file = "cons_locs_pathways_mc4_Column.h5", name = name),
         group_path) |>
    lapply(FUN = \(Adj) Matrix(Adj, sparse = TRUE)) -> Adj.mat
  cMat.info |>
    filter(group_to == groupTo[i] & !is.na(group_from)) |>
    pull(group_from) |>
    rep(times = sapply(Adj.mat, nrow)) -> rowNames
  rowNames <- paste0("NEO_", unlist(lapply(Adj.mat, \(x) seq_len(nrow(x)))), ":", rowNames)
  Adj.mat <- do.call(rbind, Adj.mat)
  rownames(Adj.mat) <- rowNames 
  colnames(Adj.mat) <- paste0("NEO_", seq_len(ncol(Adj.mat)), ":", groupTo[i])
  Adj[[i]] <- Adj.mat
}

# checking dimensions
sapply(Adj, nrow) |> unique()
sapply(Adj, ncol) |> sum()

# checking rownames
(sapply(Adj, rownames) == rownames(Adj.mat)) |> all()

# checking colnames == rownames
all(unlist(sapply(Adj, colnames)) == rownames(Adj.mat))

Adj <- do.call(cbind, Adj)
dim(Adj)

setdiff(ls(), c("Adj", "cMat.info"))

object.size(Adj) |>
  print(units = "Mb")
# 95.8 Mb it should be doable (maybe!)


##################################################################
## NETWORK COSTRUCTION

Markram_mc4 <- graph_from_adjacency_matrix(Adj, mode = "directed")
rm(Adj)

################################################################################

Markram_mc4 <- as_tbl_graph(Markram_mc4)

# check
all(
  V(Markram_mc4)$name |>
    strsplit(split = ":") |>
    sapply(\(x) x[2L]) == V.attr$Popolation
)

strsplit(V.attr$Popolation, split = "_") |>
  sapply(\(x) x[1L]) -> Layer

strsplit(V.attr$Popolation, split = "_") |>
  sapply(\(x) x[-1L]) -> Type

sapply(Type, \(x) paste(x, collapse = "_")) -> Type

# check
all(paste(Layer, Type, sep = "_") == V.attr$Popolation)
V.attr <- cbind(Layer, Type, V.attr)

V(Markram_mc4)$Layer <- V.attr$Layer
V(Markram_mc4)$Type <- V.attr$Type
V(Markram_mc4)$Popolation <- V.attr$Popolation
V(Markram_mc4)$x <- V.attr$x
V(Markram_mc4)$y <- V.attr$y
V(Markram_mc4)$z <- V.attr$z
V(Markram_mc4)$nCellAff <- V.attr$nCellAff
V(Markram_mc4)$nCellEff <- V.attr$nCellEff
V(Markram_mc4)$nSynAff <- V.attr$nSynAff
V(Markram_mc4)$nSynEff <- V.attr$nSynEff

Markram_mc4

rm(list = setdiff(ls(), "Markram_mc4"))


########################################################################

##
## EXPLORATIVE ANALYSIS AND GRAPHS
##

# nodes rapresentation


rgl::plot3d(x = V(Markram_mc4)$x, y = V(Markram_mc4)$y, z = V(Markram_mc4)$z, pch = 16, cex = 0.3, xlab = "x", ylab = "y", zlab = "z", col = Community2$membership)



# compact representation of the network

Markram_mc4 |>
  morph(to_contracted, Popolation) -> reduced_Markram_mc4
reduced_Markram_mc4 <- reduced_Markram_mc4$contracted 
reduced_Markram_mc4 |>
  activate(nodes)|>
  mutate(n.Neu = sapply(X = .N()$.tidygraph_node_index, FUN = length))|>
  activate(edges)|>
  mutate(n.edge = sapply(X = .E()$.tidygraph_edge_index, FUN = length)) -> reduced_Markram_mc4

reduced_Markram_mc4 <- delete_edge_attr(reduced_Markram_mc4, name = ".tidygraph_edge_index")
reduced_Markram_mc4 <- delete_edge_attr(reduced_Markram_mc4, name = ".orig_data")
reduced_Markram_mc4 <- delete_vertex_attr(reduced_Markram_mc4, name = ".orig_data")
reduced_Markram_mc4 <- delete_vertex_attr(reduced_Markram_mc4, name = ".tidygraph_node_index")

reduced_Markram_mc4|>
  morph(to_undirected) -> reduced_Markram_mc4.2
reduced_Markram_mc4 |>
  to_undirected() -> reduced_Markram_mc4.2

# -----
reduced_Markram_mc4
#rm(Community, Markram_mc4)
#save.image("E:/Laurea Magistrale/Tesi di laurea/Dati/reduced_Markram_mc4.RData")
# -----

plot(reduced_Markram_mc4, layout = layout.circle, vertex.label=NA)
RCy3::createNetworkFromIgraph(reduced_Markram_mc4.2)





# matrix representation
ggraph(reduced_Markram_mc4, 'matrix', sort.by=node_rank_spectral()) + 
  geom_edge_point(mirror = TRUE, show.legend = TRUE) + 
  coord_fixed()+
  theme_void()

# standard representation
ggraph(reduced_Markram_mc4, layout = 'kk', circular = TRUE) + 
  geom_edge_arc(color ="gray") + 
  coord_fixed()+
  theme_void()

ggraph(reduced_Markram_mc4, layout = 'linear') + 
  geom_edge_link() + 
  geom_node_point(size=V(reduced_Markram_mc4)$n.Neu)+
  theme_void()



# ---------- degree distributions, total, in- and out- -------------------------------

Markram_mc4 |>
  activate(nodes)|>
  mutate(In_degree = centrality_degree(mode = "in"),
         Out_degree = centrality_degree(mode = "out"),
         Degree = centrality_degree(mode = "all")) -> Markram_mc4


tibble("totDegree" = V(Markram_mc4)$Degree, 
       "inDegree" = V(Markram_mc4)$In_degree, 
       "outDegree" = V(Markram_mc4)$Out_degree, 
       "Population" = V(Markram_mc4)$Popolation) -> Degree
#pivot_longer(cols = c(inDegree, outDegree), names_to = "Type", values_to = "Degree") -> Degree


pop <- unique(Degree$Population)

for(i in pop){
  Degree |> 
    dplyr::filter(Population == i) |>
    ggplot(mapping = aes(x = Degree))+
    geom_histogram(mapping = aes(y = ..density.., fill = Type), alpha = 0.5) + 
    #geom_rug(mapping = aes(x = Degree, color = Type),lwd= 0.7, show.legend = FALSE)+
    #geom_density(mapping = aes(x = Degree, color = Type), lwd = 1)+
    theme_hc()+
    scale_color_tableau()+
    scale_fill_tableau()+
    labs(title = paste("In- and Out-degree distribution for", i,"population"))
  ggsave(filename = paste0("degree distribution_",i,".pdf"))
}

# log-log in-degree distribution
occur <- as.vector(table(V(Markram_mc4)$In_degree))
occur <- occur/sum(occur)
p <- occur/sum(occur)
y <- rev(cumsum(rev(p)))
x <- as.numeric(names(table(V(Markram_mc4)$In_degree)))
degree_log <- as_tibble('y' = y,
                        'x' = x)
degree_log |>
  ggplot(mapping = aes(x = log(x), y = log(y)))+
  geom_line(size = 1.4, color = "red")+
  theme_hc()+
  labs(title = "Log-Log scale Global In-Degree distribution")+
  ylab("log(density)")+
  xlab("log(degree)")
ggsave("log-lod inDegree.pdf")

fit_power_law(V(Markram_mc4)$In_degree) # statistical value: 0.03387035; pvalore: 0.8602261
# log-log out-degree distribution

occur <- as.vector(table(V(Markram_mc4)$Out_degree))
occur <- occur/sum(occur)
p <- occur/sum(occur)
y <- rev(cumsum(rev(p)))
x <- as.numeric(names(table(V(Markram_mc4)$Out_degree)))
degree_log <- as_tibble('y' = y,
                        'x' = x)
degree_log |>
  ggplot(mapping = aes(x = log(x), y = log(y)))+
  geom_line(size = 1.4, color = "red")+
  theme_hc()+
  labs(title = "Log-Log scale Global Out-Degree distribution")+
  ylab("log(density)")+
  xlab("log(degree)")
ggsave("log-log outdegree.pdf")

fit_power_law(V(Markram_mc4)$Out_degree) # statistical value: 0.06414975; pvalore: 0.8243975

# power law degree distribution
## 1.: goodness-of-fit test
fit_power_law(V(Markram_mc4)$Degree) # Null hypothesis rejected

# 2.: log-log degree distribution
occur <- as.vector(table(V(Markram_mc4)$Degree))
occur <- occur/sum(occur)
p <- occur/sum(occur)
y <- rev(cumsum(rev(p)))
x <- as.numeric(names(table(V(Markram_mc4)$Degree)))
degree_log <- as_tibble('y' = y,
                        'x' = x)
degree_log |>
  ggplot(mapping = aes(x = log(x), y = log(y)))+
  geom_line(size = 1.4, color = "red")+
  theme_hc()+
  labs(title = "Log-Log scale Global Degree Distribution")+
  ylab("log(density)")+
  xlab("log(degree)")
ggsave(filename = "Log-Log tot_degree distribution.pdf")

# 3.: global in/out-degree distribution
Degree |> 
  ggplot(mapping = aes(x = Degree))+
  geom_histogram(mapping = aes(y = ..density..), alpha = 0.5, color = "black") + 
  #geom_rug(mapping = aes(x = Degree, color = Type),lwd= 0.7, show.legend = FALSE)+
  #geom_density(mapping = aes(x = Degree, color = Type), lwd = 1)+
  facet_wrap(~Type)+
  theme_hc()+
  scale_color_tableau()+
  scale_fill_tableau()+
  labs(title = "Global in/out-Degree Distribution")+
  xlab("Degree")
ggsave(filename = "Global in-_out-Degree Distribution.pdf")

# 4.: global degree distribution

Degree |> 
  ggplot(mapping = aes(x = totDegree))+
  geom_histogram(mapping = aes(y = ..density..), alpha = 0.5, color = "black") + 
  #geom_rug(mapping = aes(x = Degree, color = Type),lwd= 0.7, show.legend = FALSE)+
  #geom_density(mapping = aes(x = Degree, color = Type), lwd = 1)+
  theme_hc()+
  scale_color_tableau()+
  scale_fill_tableau()+
  labs(title = "Global Degree Distribution")+
  xlab("Degree")
ggsave(filename = "Global Degree Distribution.pdf")
summary(Degree$totDegree)

# ---------- density ----------------

edge_density(Markram_mc4) # 0.0081; network is sparse



# ----------- central measures-----------------
Markram_mc4 |>
  activate(nodes) |>
  mutate(Betweennes = centrality_betweenness(),
         Clossness = centrality_closeness(),
         Eigen = centrality_eigen(directed = TRUE),
         Pagerank = centrality_pagerank()) -> Markram_mc4 # alfa = 0.85


#Problem while computing `Clossness = centrality_closeness()`. i At centrality.c:2874 :closeness centrality is not well-defined for disconnected graphs 

Centrality <- tibble(Name = V(Markram_mc4)$name,
                     Population = V(Markram_mc4)$Popolation,
                     Betweenness = V(Markram_mc4)$Betweennes,
                     Clossness = V(Markram_mc4)$Clossness,
                     Eigen = V(Markram_mc4)$Eigen,
                     Pagerank = V(Markram_mc4)$Pagerank)

# betweenness
Centrality|>
  arrange(desc(Betweenness))|>
  slice_head( n = 5)

Centrality|>
  arrange(desc(Betweenness))|>
  slice_tail(n = 5)

# closeness

Centrality|>
  arrange(desc(Clossness))|>
  slice_head( n = 5)

Centrality|>
  arrange(desc(Clossness))|>
  slice_tail( n = 5)

# eigen

Centrality|>
  arrange(desc(Eigen))|>
  slice_head( n = 5)


Centrality|>
  arrange(desc(Eigen))|>
  slice_tail( n = 5)

# pagerank

Centrality|>
  arrange(desc(Pagerank))|>
  slice_head( n = 5)

Centrality|>
  arrange(desc(Pagerank))|>
  slice_tail( n = 5)

#-------------assortativity----------------------
assortativity_degree(Markram_mc4)
assortativity(graph = Markram_mc4, types1 = V(Markram_mc4)$nCellAff, types2 =  V(Markram_mc4)$nCellEff, directed = TRUE)
# 0.1259985
#assortativity_nominal(Markram_mc4, types = V(Markram_mc4)$Type)
#assortativity_nominal(Markram_mc4, types = V(Markram_mc4)$Popolation)

#-----------------k-core decomposition-------------------------
Markram_mc4 |>
  activate(nodes) |>
  mutate(AllKCores = node_coreness(mode = "all"),
         OutKCores = node_coreness(mode = "out"),
         InKCores = node_coreness(mode = "in")) -> Markram_mc4


tibble(InKCores = V(Markram_mc4)$InKCores,
       OutKCores = V(Markram_mc4)$OutKCores)|>
  pivot_longer(cols = c(InKCores,OutKCores), names_to = "In_Out", values_to = "K_Cores")|>
  #mutate(In_Out = as.factor(In_Out))|>
  ggplot()+
  geom_histogram(mapping = aes(x = K_Cores, fill = In_Out))+
  #geom_histogram(mapping = aes(x = OutKCores),  alpha = 0.4)+
  scale_fill_tableau()+
  theme_hc()
#ggsave("k_cores_distribution.pdf")


decomposition <- coreness(Markram_mc4)



#------------- reciprocity ---------------------


reciprocity(Markram_mc4) # 0.02433

# -------- diameter ----------------------


diameter(Markram_mc4) # 7

# --------- modularity: leading_eigenvalues method ---------------- 



Markram_mc4 |>
  activate(nodes)|>
  mutate(Community = group_leading_eigen()) -> Markram_mc4

table(V(Markram_mc4)$Community)

tibble(Pop = V(Markram_mc4)$Popolation)|>
  ggplot(mapping = aes(x = Pop))+
  geom_bar(color = "gray")+
  theme(axis.text.x = element_text(angle = 70, hjust = 1))+
  theme_hc()+
  labs(title = "Frequency ditribution of neuron populations")+
  xlab("Populations")+
  ylab("Number of neurons")
ggsave(filename = "populations.pdf")

tibble(Pop = V(Markram_mc4)$Popolation,
       Comm = V(Markram_mc4)$Community)|>
  dplyr::filter(Comm == 5)|>
  ggplot(mapping = aes(x = Pop))+
  geom_bar(color = "gray")+
  theme(axis.text.x = element_text(angle = 70, hjust = 1))+
  theme_hc()+
  labs(title = "Frequency distribution of neuron populations in the 5° community")+
  xlab("Populations")+
  ylab("Number of neurons")
ggsave(filename = "5_community.pdf")

tibble(Pop = V(Markram_mc4)$Popolation,
       Comm = V(Markram_mc4)$Community)|>
  dplyr::filter(Comm == 1)|>
  ggplot(mapping = aes(x = Pop))+
  geom_bar(color = "gray")+
  theme(axis.text.x = element_text(angle = 70, hjust = 1))+
  theme_hc()+
  labs(title = "Frequency distribution of neuron populations in the 1° community")+ # la più numerosa
  xlab("Populations")+
  ylab("Number of neurons")
ggsave(filename = "1_community.pdf")


tibble(Pop = V(Markram_mc4)$Popolation,
       Comm = V(Markram_mc4)$Community)|>
  dplyr::filter(Comm == 2)|>
  ggplot(mapping = aes(x = Pop))+
  geom_bar(color = "gray")+
  theme(axis.text.x = element_text(angle = 70, hjust = 1))+
  theme_hc()+
  labs(title = "Frequency distribution of neuron populations in the 2° community")+
  xlab("Populations")+
  ylab("Number of neurons")
ggsave(filename = "2_community.pdf")

tibble(Pop = V(Markram_mc4)$Popolation,
       Comm = V(Markram_mc4)$Community)|>
  dplyr::filter(Comm == 3)|>
  ggplot(mapping = aes(x = Pop))+
  geom_bar(color = "gray")+
  theme(axis.text.x = element_text(angle = 70, hjust = 1))+
  theme_hc()+
  labs(title = "Frequency distribution of neuron populations in the 3° community")+
  xlab("Populations")+
  ylab("Number of neurons")
ggsave(filename = "3_community.pdf")

tibble(Pop = V(Markram_mc4)$Popolation,
       Comm = V(Markram_mc4)$Community)|>
  dplyr::filter(Comm == 4)|>
  ggplot(mapping = aes(x = Pop))+
  geom_bar(color = "gray")+
  theme(axis.text.x = element_text(angle = 70, hjust = 1))+
  theme_hc()+
  labs(title = "Frequency distribution of neuron populations in the 4° community")+
  xlab("populations")+
  ylab("Number of neurons")
ggsave(filename = "4_community.pdf")

rgl::plot3d(x = V(Markram_mc4)$x, y = V(Markram_mc4)$y, z = V(Markram_mc4)$z, pch = 16, cex = 0.3, xlab = "x", ylab = "y", zlab = "z", col = V(Markram_mc4)$Community)




#----------- percolation: node connectivity impact ------------
# Connectivity loss indices quantify the decrease in the number of relationships between each node when one or several components are removed. Swan_connectivty measures the loss of connectivity when excluding a node.

Markram_mc4 |>
  activate(nodes)|>
  mutate(Impact = node_connectivity_impact()) -> Markram_mc4 


#-------- adhesion and cohesion -------------

edge_connectivity(Markram_mc4) # adhesion
vertex_connectivity(Markram_mc4) # cohesion

#------- average path length --------------

mean_distance(graph = Markram_mc4) # 2.466602
?mean_distance

distance_table(Markram_mc4, directed = TRUE)

all_paths <- distances(graph = Markram_mc4)
max(all_paths)

##-----
rm(Community)

#------- ERG Modelling---------------------

subgraph5 <- subgraph(Markram_mc4, V(Markram_mc4)$Community ==5)
coords <- tibble("x" = V(subgraph5)$x,
                 "y" = V(subgraph5)$y,
                 "z" = V(subgraph5)$z)
#n <- length(V(subgraph5)$x)

dist(coords, diag = TRUE, upper = TRUE) |> as.matrix() -> Dst
#V(subgraph5)$Distance <- Dst
E(subgraph5)$Distance <- Dst
V(subgraph5)$Idegree <- subgraph5 |> degree(mode = "in")
V(subgraph5)$Odegree <- subgraph5 |> degree(mode = "out")


subgraph5 <- asNetwork(subgraph5) 
set.network.attribute(subgraph5, attrname = "Distance", value = Dst)
rm(Markram_mc4, coords, Dst)



get.network.attribute(x = subgraph5, attrname = "Distance") -> Dist


# MODELS--------------------

# Initial estimate of terms

subgraph5.Summary <- summary(subgraph5 ~ mutual + asymmetric + m2star + isolates + gwidegree + gwodegree + ttriple + ctriple)

subgraph5.Summary



#----- Model Fitting ------


# Bernoulli independence: exponential model
system.time(modello0 <- ergm(subgraph5 ~ edges , estimate = "MPLE"))
summary(modello0) # 712 mb
# AIC: 3575625  BIC: 3575639

# Bernoulli independence + distanza
modello1 <- ergm(subgraph5 ~ edges + dyadcov("Distance"), estimate = "MPLE")
summary(modello1)
# AIC: 3428363  BIC: 3428420

# Dyad independence + mutual
system.time(modello2 <- ergm(subgraph5 ~ edges  + mutual , estimate = "MPLE"))
summary(modello2)
# AIC: 3427914  BIC: 3427986


#Model 2.2
modello2 <- ergm(subgraph5 ~ edges  + mutual + gwidegree + gwodegree, estimate = "MPLE")
summary(modello2)
# AIC: 3645371  BIC: 3645457


# Markov models
#modello3 <- ergm(subgraph5 ~ edges  + mutual  + dyadcov("Distance") + istar(2) + ostar(2), estimate = "MPLE")
system.time(modello3 <- ergm(subgraph5 ~ edges + mutual  + istar(2) + ostar(2), estimate = "MPLE"))
summary(modello3)
# AIC: 3356053  BIC: 3356111

# Markov models: model 4
modello4 <- ergm(subgraph5 ~ edges  + mutual  + m2star + istar(2) + ostar(2), estimate = "MPLE")
summary(modello4)


# Markov models: model 5
modello5 <- ergm(subgraph5 ~ edges  + mutual + istar(2) + ostar(2)+ gwidegree + gwodegree, estimate = "MPLE")
summary(modello5)


# Markov model: model 6
system.time(modello6 <- ergm(subgraph5~ edges + istar(2) + ostar(2) + istar(3) + ostar(3) + istar(4) + ostar(4) + mutual, estimate = "MPLE"))
summary(modello6)







# AICc calculation
n <- 399342

# (2*n*k)/(n-K-1)

# modello0
3575625 - 2 + (2 * n * 1)/(n - 1-1)

# modello1 
3431173 - 4 + (2*2*n)/(n - 2 - 1)
# modello2
3481760 - 4 + (2*2*n)/(n-2-1)

# modello3
3413172 - 6 + (2*3*n)/(n - 3 - 1)

# modello4
3571208 - 2 + (2*2*n)/(n-2-1)

# modello5
3430539 - 6 + (2*n*3)/(n-3-1)

# modello6
3356053 - 8 + (2*n*4)/(n-4-1)

# modello7
3332040 - 16 + (2*n*8)/(n-8-1)





#------------ Evaluation of Goodness-of-Fit -----------------
subgraph5 <- asIgraph(subgraph5)


sim_100 <- simulate(modello3, nsim = 100, output = "edgelist", verbose = 2, control = control.simulate.formula(MCMC.burnin = 1000, MCMC.interval = 100)) 
sim_100 <- lapply(sim_100, FUN = asIgraph)
gc()


In <- matrix(data = NA, nrow = 3623, ncol = 100)
Out <- matrix(data = NA, nrow = 3623, ncol = 100)


for(i in 1:100){
  In[,i] <- degree(sim_100[[i]], mode = "in")
  Out[,i] <- degree(sim_100[[i]], mode = "out")
}

# In-Degree
apply(In, MARGIN = 2, FUN = function(x) as.vector(prop.table(table(x)))[1:200]) -> In
In <- t(In)
IN <- apply(In, MARGIN = 2, FUN = median)


pdf(file = "in-degree-gof.pdf")
boxplot(In, xlab = "in-degree", ylab = "proportion of nodes")
#points(as.vector(prop.table(table(degree(subgraph5, mode = "in"))))[1:200], col = 2, cex = 1, pch=".")
lines(as.vector(prop.table(table(degree(subgraph5, mode = "in"))))[1:200], col = 2, lwd = 2)
dev.off()


# Out-Degree
apply(Out, MARGIN = 2, FUN = function(x) as.vector(prop.table(table(x)))[1:200]) -> Out
Out <- t(Out)
pdf(file = "out-degree-gof.pdf")
boxplot(Out, xlab = "out-degree", ylab = "proportion of nodes")
#points(as.vector(prop.table(table(degree(subgraph5, mode = "out"))))[1:200], col = 2, cex = 1, pch=".")
lines(as.vector(prop.table(table(degree(subgraph5, mode = "out"))))[1:200], col = 2, lwd = 2)
dev.off()



#### alternative graphs for the GOF #######################ù

#----------- In- and Out degree distributions of simulated networks

In <- matrix(data = NA, nrow = 3623, ncol = 100)
Out <- matrix(data = NA, nrow = 3623, ncol = 100)


for(i in 1:100){
  In[,i] <- degree(sim_100[[i]], mode = "in")
  Out[,i] <- degree(sim_100[[i]], mode = "out")
}


apply(In, MARGIN = 2, FUN = function(x) as.vector(prop.table(table(x)))[1:3623]) -> In
apply(In, MARGIN = 2, FUN = median) -> In



apply(Out, MARGIN = 2, FUN = function(x) as.vector(prop.table(table(x)))) -> Out
Out <- t(Out)
apply(Out, MARGIN = 2, FUN = median) -> Out

dt <- cbind("In" = In, "OUt" = Out, "In_Obs" = as.vector(prop.table(table(degree(subgraph5, mode = "in")))), "Out_Obs" = as.vector(prop.table(table(degree(subgraph5, mode = "out")))))
dt <- as.data.frame(dt)
dt <- 
  
  
  dt |>
  ggplot()|>
  geom_line



#--------- Geodesic Distance of simulated networks


GeoD <- matrix(0, nrow = 6, ncol = 100)
GeoNA <- rep(NA, 100)

for(i in 1:100){
  d <- distance_table(sim_100[[i]])
  GeoNA[i] <- d$unconnected
  if(length(d$res) < 6){
    dif <- 6 - length(d$res)
    GeoD[,i] <- c(d$res, rep(0, dif))
  }
  else{
    GeoD[,i] <- d$res
  }
}

GeoDist <- rbind(GeoD,  GeoNA)

GeoDist2 <- apply(X = GeoDist, MARGIN = 1, FUN = function(x) x/sum(x))
colnames(GeoDist2) <- c(1,2,3,4,5,6,"NA")
boxplot(GeoDist2, xlab="minimum geodesic distance", ylab="proportion of dyads", ylim=c(0, 0.80))

observed_geod <- distance_table(subgraph5)
observed_geod <- c(observed_geod$res, observed_geod$unconnected)
observed_geod <- observed_geod/sum(observed_geod)
lines(observed_geod, col = 2, lwd = 2)
