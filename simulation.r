library(igraph)
library(ggplot2)
library(colorspace)

# Stochastic block model:
sim_graph <- sample_sbm(200, pref.matrix = pm, block.sizes = rep(20, times = 10))


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

# Our end frame will have the patient ID, whether they are infected and the iteration
# on which they are infected.
# Infection should radiate outward from the idnex.
sim_patients <- patient_frame |>
    mutate(infected = FALSE, iteration = NA)

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
    
    # This gives me the connected nodes from the target nodes. Note the distance
    # is always one because we want the next iteration to be connected to 
    # the previous
    possible_infected <- distances(sim_graph, v = patients_spreading_from) == 1
    possible_infected <- colnames(possible_infected)[colSums(possible_infected) > 0]

    # Infection probability decays exponentially with the number of hops.
    are_infected <- runif(length(possible_infected)) < baseline_chance^i

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
base_detection_chance <- 0.6

sim_patients <- sim_patients |>
    mutate(detected = NA, iteration_detected = NA) |>
    mutate(
        detected = if_else(sID == index_patient, TRUE, FALSE),
        iteration_detected = if_else(sID == index_patient, 0, NA)
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
    are_infected <- runif(length(possible_detected)) < base_detection_chance^i

    sim_patients <- sim_patients |>
        mutate(
            iteration_detected = if_else(sID %in% possible_detected[are_infected] & !detected, i, iteration),
            detected = if_else(sID %in% possible_detected[are_infected] & !detected, TRUE, detected)
        )
}

overall_detected <- sim_patients |>
    filter(detected) |>
    pull(sID)

sim_graph <- sim_graph |>
    set_vertex_attr("detected", value = V(sim_graph)$name %in% overall_detected) |>
    set_vertex_attr("iteration_detected", value = sim_patients$iteration_detected, index = sim_patients$sID)

sim_graph |>
    #    delete_vertices(V(sim_graph)$infected == FALSE) |>
    ggplot(aes(x = x, y = y, xend = xend, yend = yend)) +
    geom_edges() +
    geom_nodes(aes(colour = factor(iteration), shape = detected), size = 4) +
    scale_colour_discrete_sequential(rev = FALSE, palette = "reds", na.value = "lightgrey") +
    labs(colour = "Iteration of infection")
