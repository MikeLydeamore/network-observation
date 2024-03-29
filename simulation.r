library(igraph)
library(ggplot2)
library(colorspace)
library(tidyverse)
library(ggraph)
library(intergraph)
library(ggnetwork)


set.seed(462612)
pm <- matrix(0.00002, nrow = 45, ncol = 45)
diag(pm) <- 0.4

set.seed(462612)

# Stochastic block model:
sim_graph <- sample_sbm(4500, pref.matrix = pm, block.sizes = rep(100, times = 45))

sim_graph <- sim_graph |>
    set_vertex_attr("name", value = paste0("s", 1:4500))


# Clusters will act as our proxy for wards. Patients are more likely to be connected
# within a cluster than between so this should suffice.
clusters <- cluster_infomap(sim_graph)

sim_graph <- sim_graph |>
    set_vertex_attr("cluster", value = as.numeric(membership(clusters)))

# We will pick an index ward and then a patient within it:
ward <- sample(unique(vertex_attr(sim_graph, "cluster")), size = 1)

# Convenient structure that has the patient and which ward they are in.
patient_frame <- tibble(
    sID = V(sim_graph)$name,
    ward = membership(clusters)
)

index_patient <- patient_frame |>
    filter(ward == ward) |>
    sample_n(1) |>
    pull(sID)

index_ward <- patient_frame |>
    filter(sID == index_patient) |>
    pull(ward)

# Our end frame will have the patient ID, whether they are infected and the iteration
# on which they are infected.
# Infection should radiate outward from the idnex.
sim_patients <- patient_frame |>
    mutate(infected = FALSE, iteration = NA)

initial_num_infected <- sim_patients |>
  filter(ward == index_ward) |>
  nrow()

sim_patients <- sim_patients |>
    mutate(
        infected = ifelse(ward == index_ward, runif(n = initial_num_infected) < 0.1, FALSE),
        iteration = if_else(infected, 0, NA)
    )

iterations <- 5
baseline_chance <- 0.5
for (i in 1:iterations) {
    patients_spreading_from <- sim_patients |>
        filter(infected, iteration == (i - 1)) |>
        pull(sID)
    
    # This gives me the connected nodes from the target nodes. Note the distance
    # is always one because we want the next iteration to be connected to 
    # the previous
    possible_infected <- distances(sim_graph, v = patients_spreading_from) == 1
    possible_infected <- colnames(possible_infected)[colSums(possible_infected) > 0]

    # Infection probability decays exponentially with the number of hops.
    are_infected <- runif(length(possible_infected)) < baseline_chance

    sim_patients <- sim_patients |>
        mutate(
            iteration = if_else(sID %in% possible_infected[are_infected] & !infected, i, iteration),
            infected = if_else(sID %in% possible_infected[are_infected] & !infected, TRUE, infected)
        )
}

# It's really hard to get the attributes onto the graph so IMO it's just easier
# to do it at the end.
overall_infected <- sim_patients |>
    filter(infected) |>
    pull(sID)
sim_graph <- sim_graph |>
    set_vertex_attr("infected", value = V(sim_graph)$name %in% overall_infected) |>
    set_vertex_attr("iteration", value = sim_patients$iteration, index = sim_patients$sID)

sim_graph |>
    #    delete_vertices(V(sim_graph)$infected == FALSE) |>
    ggplot(aes(x = x, y = y, xend = xend, yend = yend)) +
    geom_edges() +
    geom_nodes(aes(colour = factor(iteration)), size = 4) +
    scale_colour_discrete_sequential(rev = FALSE, palette = "reds", na.value = "lightgrey") +
    labs(colour = "Iteration of infection")

## Detection:

iterations <- iterations

base_detection_chance <- 0.8

sim_patients <- sim_patients |>
    mutate(detected = NA, iteration_detected = NA) |>
    mutate(
        detected = if_else(ward == index_ward, TRUE, FALSE),
        iteration_detected = if_else(ward == index_ward, 0, NA)
    )

for (i in 1:iterations) {
    patients_spreading_from <- sim_patients |>
        filter(infected, iteration == (i - 1)) |>
        pull(sID)

    # This gives me the connected nodes from the target nodes. Note the distance
    # is always one because we want the next iteration to be connected to
    # the previous
    possible_detected <- distances(sim_graph, v = patients_spreading_from) == 1
    possible_detected <- colnames(possible_detected)[colSums(possible_detected) > 0]

    # Infection probability decays exponentially with the number of hops.
    are_detected <- runif(length(possible_detected)) < base_detection_chance^i

    sim_patients <- sim_patients |>
        mutate(
            iteration_detected = if_else(sID %in% possible_detected[are_detected] & !detected, i, iteration_detected),
            detected = if_else(sID %in% possible_detected[are_detected] & !detected, TRUE, detected)
        )
    sim_patients |> filter(detected, is.na(iteration_detected))
}

overall_detected <- sim_patients |>
    filter(detected) |>
    pull(sID)

sim_graph <- sim_graph |>
    set_vertex_attr("detected", value = V(sim_graph)$name %in% overall_detected) |>
    set_vertex_attr("iteration_detected", value = sim_patients$iteration_detected, index = sim_patients$sID) |>
    set_vertex_attr("tested", value = !is.na(sim_patients$iteration_detected)) |>
    set_vertex_attr("index_patient", value = V(sim_graph)$name == index_patient)

# sim_graph |>
#     #    delete_vertices(V(sim_graph)$infected == FALSE) |>
#     ggplot(aes(x = x, y = y, xend = xend, yend = yend)) +
#     geom_edges() +
#     geom_nodes(aes(colour = factor(iteration), shape = detected), size = 4) +
#     scale_colour_discrete_sequential(rev = FALSE, palette = "reds", na.value = "lightgrey") +
#     labs(colour = "Iteration of infection")
# 
# 
# ggnetwork(sim_graph, by = "infected") |>
#     ggplot(aes(x = x, y = y, xend = xend, yend = yend)) +
#         geom_edges(colour = "grey") +
#         geom_nodes(aes(fill = factor(cluster), shape = infected), size = 4, alpha = 0.6) +
#         geom_nodes(aes(colour = tested), size = 4, shape = 4, stroke = 1) +
#         # scale_colour_discrete_sequential(rev = FALSE, palette = "reds", na.value = "lightgrey") +
#         labs(colour = "Tested", fill = "Ward") +
#         scale_colour_manual(values = c("black", "red")) +
#         scale_shape_manual(values = c(21, 24)) +
#         facet_wrap(~infected)
# 
# 
# sim_graph |>
#     delete_vertices(V(sim_graph)$infected == FALSE) |>
#     ggplot(aes(x = x, y = y, xend = xend, yend = yend)) +
#     geom_edges(colour = "grey") +
#     geom_nodes(data = function(x) {
#         x[x$index_patient, ]
#     }, size = 10, colour = "red", alpha = 0.4) +
#     geom_nodes(aes(fill = factor(cluster), shape = infected), size = 4, alpha = 0.6) +
#     geom_nodes(aes(colour = tested), size = 4, shape = 4, stroke = 1) +
#     # scale_colour_discrete_sequential(rev = FALSE, palette = "reds", na.value = "lightgrey") +
#     labs(colour = "Tested", fill = "Ward") +
#     scale_colour_manual(values = c("black", "red")) +
#     scale_shape_manual(values = c(21, 24)) +
#     theme_void()


# sim_layout <- get_layout_radial(sim_patients, sim_graph)
    
# 
# ggnetwork(sim_graph, layout = as.matrix(sim_layout[, c("x", "y")])) |>
#     ggplot(aes(x = x, y = y, xend = xend, yend = yend)) +
#     #geom_edges(colour = "lightgrey") +
#     geom_nodes(aes(shape = infected), size = 4, alpha = 0.6, fill = "black") +
#     geom_nodes(aes(colour = tested), size = 4, shape = 4, stroke = 1) +
#     labs(colour = "Tested", fill = "Ward") +
#     scale_colour_manual(values = c("black", "red")) +
#     scale_shape_manual(values = c(21, 24)) +
#     coord_fixed() +
#     theme_void()

ward_network <- contract(
    sim_graph,
    mapping = membership(clusters),
    vertex.attr.comb = function(x) {
        patient_frame |>
            filter(sID %in% x) |>
            pull(ward) |>
            unique()
    }
)

vertex_attr(ward_network, "name")

ward_adj_strength <- as_adj(ward_network)
ward_adj_strength <- as.matrix(ward_adj_strength)
ward_adj <- ward_adj_strength
diag(ward_adj) <- 0

ward_adj[ward_adj > 0] <- 1

ward_adj <- Matrix::drop0(ward_adj)

is_connected(ward_network)


