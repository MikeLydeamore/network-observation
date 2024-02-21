library(HospitalNetwork)
library(lubridate)
library(ggraph)

# Generate a synthetic network
# sID = Patient
# fID = Facility
# Adate = Admission date
# Dddate = Discharge date
n_subjects <- 1000000
synthetic_db <- HospitalNetwork::create_fake_subjectDB_clustered(
    n_subjects = n_subjects,
    n_facilities = 300,
    avg_n_stays = 10,
    n_clusters = 10
)

# Patients who have an infection:

num_infected_patients <- round(n_subjects * 0.00005)
infected_patients <- sample(unique(synthetic_db$sID), num_infected_patients)

dates_of_infection <- synthetic_db %>%
    filter(sID %in% infected_patients) %>%
    group_by(sID) %>%
    arrange(sID, Adate) %>%
    slice_sample() %>%
    rowwise() %>%
    mutate(
        length_of_stay = as.numeric(difftime(Ddate, Adate, units = "days")),
        date_of_infection = Adate + days(sample(0:length_of_stay, 1))
    ) %>%
    select(sID, date_of_infection)

perfect_detection <- synthetic_db %>%
    left_join(dates_of_infection, by="sID") %>%
    mutate(
        infected = !is.na(date_of_infection) & Adate >= date_of_infection,
        detected = infected
    )

uniform_detection <- synthetic_db %>%
    left_join(dates_of_infection, by="sID") %>%
    mutate(
        infected = !is.na(date_of_infection) & Adate >= date_of_infection
    ) %>%
    group_by(sID) %>%
    mutate(
        detected_admission_number = 
        ifelse(
            length(which(infected)) == 0, NA,
            ifelse(which(infected)[1] == n(), n(),
                sample(which(infected)[1]:n(), size=1)
            )   
        )
    ) %>%
    mutate(detected = row_number() >= detected_admission_number) %>%
    select(-detected_admission_number)

write.csv(perfect_detection, "data/perfect_detection.csv", row.names=F)
write.csv(uniform_detection, "data/uniform_detection.csv", row.names = F)

detection_probs <- tibble(
    fID = synthetic_db$fID %>%
        unique()
) |>
    mutate(
        detection_prob = c(rep(1, 10), rep(0.1, n()-10))
    )
targetted_detection <- synthetic_db %>%
    left_join(dates_of_infection, by = "sID") %>%
    mutate(infected = !is.na(date_of_infection) & Adate >= date_of_infection) %>%
    left_join(detection_probs, by = "fID") %>%
    mutate(
        detected = if_else(
            infected,
            runif(n()) < detection_prob,
            FALSE
        )
    ) %>%
    group_by(sID) %>%
    mutate(
        first_detected = if_else(
            any(detected),
            which(detected)[1],
            NA
        )
    ) |>
    mutate(
        new_detected = row_number() >= first_detected
    )

targetted_detection |>
    select(-detected, -first_detected) |>
    rename(detected = new_detected) |>
    write.csv("data/targetted_detection.csv", row.names = F)

# I'll just get the clusters using HospitalNetwork so it's consistent.
library(HospitalNetwork)

synthetic_base <- synthetic_db |>
    checkBase() |>
    hospinet_from_subject_database()

clustered_detection_probs <- synthetic_base$cluster_infomap |>
    mutate(detection_prob = if_else(
        cluster_infomap == 4,
        1,
        0.1
    )) |>
    rename(fID = node)

clustered_detection <- synthetic_db %>%
    left_join(dates_of_infection, by = "sID") %>%
    mutate(infected = !is.na(date_of_infection) & Adate >= date_of_infection) %>%
    left_join(clustered_detection_probs, by = "fID") %>%
    mutate(
        detected = if_else(
            infected,
            runif(n()) < detection_prob,
            FALSE
        )
    ) %>%
    group_by(sID) %>%
    mutate(
        first_detected = if_else(
            any(detected),
            which(detected)[1],
            NA
        )
    ) |>
    mutate(
        new_detected = row_number() >= first_detected
    )

clustered_detection |>
    select(-detected, -first_detected) |>
    rename(detected = new_detected) |>
    write.csv("data/clustered_detection.csv", row.names = F)


clustered_detection_probs <- synthetic_base$cluster_infomap |>
    mutate(detection_prob = if_else(
        cluster_infomap == 4,
        1,
        0.0
    )) |>
    rename(fID = node)

clustered_detection <- synthetic_db %>%
    left_join(dates_of_infection, by = "sID") %>%
    mutate(infected = !is.na(date_of_infection) & Adate >= date_of_infection) %>%
    left_join(clustered_detection_probs, by = "fID") %>%
    mutate(
        detected = if_else(
            infected,
            runif(n()) < detection_prob,
            FALSE
        )
    )

clustered_detection |>
    write.csv("data/clustered_detection_2.csv", row.names = F)


normalised <- t(apply(synthetic_base$matrix, 1, function(x) {
    x / sum(x)
}))

seed_node <- "f04"

infected_nodes <- sample(rownames(normalised), size = 50, prob = normalised[seed_node, ])


synthetic_base$igraph |>
    set_vertex_attr("infected", value = V(synthetic_base$igraph)$name %in% infected_nodes) |>
    write_rds( "data/network_spread_graph.rds")


synthetic_db <- HospitalNetwork::create_fake_subjectDB_clustered(
    n_subjects = 100,
    n_facilities = 10,
    avg_n_stays = 3,
    n_clusters = 10
)


patient_edgelist <- lapply(unique(synthetic_db$fID), function(facility) {
    facility_data <- filtered |>
        filter(fID == facility)

    subjects <- unique(facility_data$sID)
    edges <- expand.grid(subjects, subjects) |>
        mutate(fID = facility)

    return(edges)
}) |>
    bind_rows() |>
    as_tibble() |>
    unique()

patient_graph <- graph_from_edgelist(as.matrix(patient_edgelist[, 1:2]), directed = FALSE)

patient_graph <- patient_graph |>
    set_edge_attr("ward", value = patient_edgelist$fID)
# Draw as multi-level network:

multilevel <- synthetic_db |>
    select(sID, fID) |>
    as.matrix() |>
    graph_from_edgelist(directed = FALSE)

multilevel <- multilevel |>
    set_vertex_attr("type",
        value = if_else(
            V(multilevel)$name %in% synthetic_db$sID,
            "Patient",
            "Ward"
        )
    )

ggplot(multilevel, aes(x = x, y = y, xend = xend, yend = yend)) +
    geom_edges(alpha = 0.3) +
    geom_nodes(aes(colour = type, shape = type))


## Pick a ward, and then a patient within it:

ward <- sample(synthetic_db$fID, size = 1)
index_patient <- synthetic_db |>
    filter(fID == ward) |>
    reframe(sample(sID, size = 1)) |>
    pull()

baseline_chance <- 0.8
iterations <- 20

infected_patients <- index_patient
adj <- as_adjacency_matrix(patient_graph)

max(distances(patient_graph))

