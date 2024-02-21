library(ergm)
library(tidyverse)
library(HospitalNetwork)
library(igraph)
library(ggnetwork)
library(intergraph)

perfect_detection <- read_csv("data/perfect_detection.csv")

add_infection_status <- function(.data, infection_status, facility) {
    .data |>
        group_by({{ facility }}) |>
        summarise(experienced_infection = any({{ infection_status }}))
}

facilities_experienced_infections <- perfect_detection |>
    add_infection_status(detected, fID)

network <- perfect_detection |>
    checkBase() |>
    hospinet_from_subject_database()

add_infection_status_to_igraph <- function(igraph, infected_facilities, experienced_infection, facility) {
    infected_facilities <- infected_facilities |>
        filter({{ experienced_infection }}) |>
        pull({{ facility }})
    
    nodes_experienced_infection = names(V(igraph)) %in%
     infected_facilities
    
    igraph |>
        set_vertex_attr(
            name = "infected", 
            value = nodes_experienced_infection
        )
}

add_clusters_to_igraph <- function(igraph, cluster_infomap) {
    igraph |>
        set_vertex_attr(
            name = "cluster", 
            value = as.numeric(cluster_infomap$cluster_infomap),
            index = cluster_infomap$node
        )
}

infection_graph <- network$igraph |>
    add_infection_status_to_igraph(facilities_experienced_infections, experienced_infection, fID) |>
    add_clusters_to_igraph(network$cluster_infomap)

ggplot(infection_graph, aes(x = x, y = y, xend = xend, yend = yend)) +
    geom_edges(alpha = 0.01) +
    geom_nodes(aes(colour = infected))

network_intergraph <- asNetwork(infection_graph)

# Fit a network accounting for the infection.
ergm_1 <- ergm(network_intergraph ~ edges + nodefactor("infected") + nodefactor("cluster"))

summary(ergm_1)

uniform_detection <- read_csv("data/uniform_detection.csv")

uniform_detection_facilities <- uniform_detection |>
    add_infection_status(detected, fID)

uniform_network <- uniform_detection |>
    checkBase() |>
    hospinet_from_subject_database()

uniform_infection_graph <- uniform_network$igraph |>
    add_infection_status_to_igraph(uniform_detection_facilities, experienced_infection, fID) |>
    add_clusters_to_igraph(uniform_network$cluster_infomap)

uniform_intergraph <- asNetwork(uniform_infection_graph)

uniform_ergm <- ergm(uniform_intergraph ~ edges + nodefactor("infected") + nodefactor("cluster"))

summary(uniform_ergm)

targetted_detection <- read_csv("data/targetted_detection.csv")

targetted_detection_facilities <- targetted_detection |>
    add_infection_status(detected, fID)

targetted_network <- targetted_detection |>
    checkBase() |>
    hospinet_from_subject_database()

targetted_infection_graph <- targetted_network$igraph |>
    add_infection_status_to_igraph(targetted_detection_facilities, experienced_infection, fID) |>
    add_clusters_to_igraph(targetted_network$cluster_infomap)

targetted_intergraph <- asNetwork(targetted_infection_graph)

targetted_ergm <- ergm(targetted_intergraph ~ edges + nodefactor("infected") + nodefactor("cluster"))

summary(targetted_ergm)

clustered_detection <- read_csv("data/clustered_detection.csv")

clustered_detection_facilities <- clustered_detection |>
    add_infection_status(detected, fID)

clustered_network <- clustered_detection |>
    checkBase() |>
    hospinet_from_subject_database()

clustered_infection_graph <- clustered_network$igraph |>
    add_infection_status_to_igraph(clustered_detection_facilities, experienced_infection, fID) |>
    add_clusters_to_igraph(clustered_network$cluster_infomap)

clustered_intergraph <- asNetwork(clustered_infection_graph)

clustered_ergm <- ergm(clustered_intergraph ~ edges + nodefactor("infected") + nodefactor("cluster"))

summary(clustered_ergm)

ggplot(clustered_infection_graph, aes(x = x, y = y, xend = xend, yend = yend)) +
    geom_edges(alpha = 0.01) +
    geom_nodes(aes(colour = infected))



clustered_detection_2 <- read_csv("data/clustered_detection_2.csv")

clustered_detection_2_facilities <- clustered_detection_2 |>
    add_infection_status(detected, fID)

clustered_network_2 <- clustered_detection_2 |>
    checkBase() |>
    hospinet_from_subject_database()

clustered_infection_graph_2 <- clustered_network_2$igraph |>
    add_infection_status_to_igraph(clustered_detection_2_facilities, experienced_infection, fID) |>
    add_clusters_to_igraph(clustered_network_2$cluster_infomap)

clustered_intergraph_2 <- asNetwork(clustered_infection_graph_2)

clustered_ergm_2 <- ergm(clustered_intergraph_2 ~ edges + nodefactor("infected") + edgecov("weight"))

summary(clustered_ergm_2)

ggplot(clustered_infection_graph_2, aes(x = x, y = y, xend = xend, yend = yend)) +
    geom_edges(alpha = 0.01) +
    geom_nodes(aes(colour = infected))


network_spread_graph <- readRDS("data/network_spread_graph.rds")

clusters <- cluster_infomap(network_spread_graph)

network_spread_graph <- network_spread_graph |>
    set_vertex_attr("cluster", value = as.numeric(membership(clusters)))

network_spread_network <- asNetwork(network_spread_graph)

network_spread_ergm <- ergm(network_spread_network ~ edges + nodefactor("infected"))
summary(network_spread_ergm)

network_spread_ergm_cluster <- ergm(network_spread_network ~ edges + nodefactor("infected") + nodefactor("cluster"))
summary(network_spread_ergm_cluster)

ggplot(network_spread_graph, aes(x = x, y = y, xend = xend, yend = yend)) +
    geom_edges(alpha = 0.01) +
    geom_nodes(aes(colour = infected))

pm <- matrix(0.002, nrow = 10, ncol = 10)
diag(pm) <- 0.4

sim_graph <- sample_sbm(200, pref.matrix = pm, block.sizes = rep(20, times = 10))

clusters <- cluster_infomap(sim_graph)

sim_graph <- sim_graph |>
    set_vertex_attr("cluster", value = as.numeric(membership(clusters)))

ggplot(sim_graph, aes(x = x, y = y, xend = xend, yend = yend)) +
    geom_edges() +
    geom_nodes(aes(colour = factor(cluster)))

sim_network <- asNetwork(sim_graph)

sim_ergm <- ergm(sim_network ~ edges + nodefactor("cluster"))

summary(sim_ergm)

get_infected_nodes <- function(graph, n_infected_nodes, seed_node = NULL) {
    adj_matrix <- as_adjacency_matrix(graph)

    normalised <- t(apply(adj_matrix, 1, function(x) {
        x / sum(x)
    }))


    if (is.null(seed_node)) {
        seed_node <- sample(rownames(normalised),  size = 1)
    }
    
    infected_nodes <- sample(rownames(normalised), size = n_infected_nodes, prob = normalised[seed_node, ])

    return(infected_nodes)
}

V(sim_graph)$name <- paste0("f", 1:200)
infected_nodes <- get_infected_nodes(sim_graph, 7, "f135")
sim_graph <- sim_graph |>
        set_vertex_attr("infected", value = V(sim_graph)$name %in% infected_nodes)

sim_ergm2 <- ergm(asNetwork(sim_graph)~edges + nodefactor("infected") + nodefactor("cluster"))

summary(sim_ergm2)

ggplot(sim_graph, aes(x = x, y = y, xend = xend, yend = yend)) +
    geom_edges() +
    geom_nodes(aes(colour = infected))


ward <- sample(unique(vertex_attr(sim_graph, "cluster")), size = 1)

patient_frame <- tibble(
    sID = V(sim_graph)$name,
    ward = membership(clusters)
)

index_patient <- patient_frame |>
    filter(ward == ward) |>
    sample_n(1) |>
    pull(sID)

distance <- distances(sim_graph)

sim_patients <- patient_frame |>
    mutate(infected = FALSE, iteration_infected = NA)

sim_patients <- sim_patients |>
    mutate(
        infected = if_else(sID == index_patient, TRUE, FALSE),
        iteration = if_else(sID == index_patient, 0, NA)
    )

iterations <- 5
baseline_chance <- 0.8
for (i in 1:iterations) {
    patients_spreading_from <- sim_patients |>
        filter(infected, iteration == (i - 1)) |>
        pull(sID)
    possible_infected <- distances(sim_graph, v = patients_spreading_from) == 1
    possible_infected <- colnames(possible_infected)[colSums(possible_infected) > 0]

    are_infected <- runif(length(possible_infected)) < baseline_chance^i

    sim_patients <- sim_patients |>
        mutate(
            infected = if_else(sID %in% possible_infected[are_infected], TRUE, infected),
            iteration = if_else(sID %in% possible_infected[are_infected], i, iteration)
        )
}

overall_infected <- sim_patients |>
    filter(infected) |>
    pull(sID)
sim_graph <- sim_graph |>
    set_vertex_attr("infected", value = V(sim_graph)$name %in% overall_infected)

sim_graph |>
#    delete_vertices(V(sim_graph)$infected == FALSE) |>
    ggplot(aes(x = x, y = y, xend = xend, yend = yend)) +
        geom_edges() +
        geom_nodes(aes(colour = infected))
