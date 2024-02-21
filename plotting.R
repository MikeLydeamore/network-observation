get_layout_radial <- function(sim_patients, sim_graph) {
    distance_matrix <- distances(sim_graph)

    index_patient <- sim_patients |>
        filter(iteration == 0) |>
        pull(sID)

    max_distance <- max(distance_matrix)

    distance_from_index <- tibble(
        sID = colnames(distance_matrix),
        distance = distance_matrix[index_patient, ]
    )

    distance_from_index |>
        group_by(distance) |>
        mutate(
            x = distance * sin(2 * row_number() * pi / n()),
            y = distance * cos(2 * row_number() * pi / n())
        )
    
}
