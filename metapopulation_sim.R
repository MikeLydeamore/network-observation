library(readr)
library(dplyr)
library(HospitalNetwork)
library(tidyr)

patient_data <- read_csv("data/perfect_detection.csv")

hospinet <- checkBase(patient_data) %>% hospinet_from_subject_database()
network <- hospinet$igraph

dt <- 1
max_t <- 365

start_day <- 0
calendar_start <- min(patient_data$Adate)

patient_data <- patient_data %>%
    mutate(
        admission_day = as.numeric(difftime(Adate, calendar_start, units = "days")),
        discharge_day = as.numeric(difftime(Ddate, calendar_start, units = "days"))
    ) %>%
    arrange(sID, Adate) %>%
    group_by(sID) %>%
    mutate(
        is_transfer = c(lag(Adate) == Adate)
    )


facilities <- patient_data %>%
    pull(fID) %>%
    unique()

metapop_data <- lapply(seq(start_day, max_t, by = dt), function(current_day) {
    print(current_day)

    colonised_t <- patient_data %>%
        filter(admission_day <= current_day & discharge_day > current_day) %>%
        filter(infected) %>%
        group_by(fID) %>%
        summarise(colonised_t = length(unique(sID)))

    admissions_t <- patient_data %>%
        filter(admission_day == current_day) %>%
        filter(!is_transfer) %>%
        group_by(fID) %>%
        summarise(admissions_t = length(unique(sID)))

    discharges_t <- patient_data %>%
        filter(discharge_day == current_day) %>%
        filter(!is_transfer) %>%
        group_by(fID) %>%
        summarise(discharges_t = length(unique(sID)))

    capacity_t <- patient_data %>%
        filter(admission_day <= current_day & discharge_day > current_day) %>%
        group_by(fID) %>%
        summarise(capacity_t = length(unique(sID)))

    tibble(fID = facilities) %>%
        full_join(colonised_t, by = "fID") %>%
        full_join(admissions_t, by = "fID") %>%
        full_join(discharges_t, by = "fID") %>%
        full_join(capacity_t, by = "fID") %>%
        mutate(day = current_day)
}) %>%
    bind_rows() %>%
    mutate(
        across(everything(), ~ replace_na(.x, 0))
    ) %>%
    arrange(day, fID)

sim_t <- start_day + dt
beta <- 1.3
gamma <- 0.01
delta <- 1

initial_infections <- paste0("f", sample(1:100, size = 20))

colonised_t <- tibble(fID = facilities, model_colonised_t = 0) %>%
    mutate(
        model_colonised_t = ifelse(fID %in% initial_infections, 1, 0)
    )

sim_results <- list()
sim_t <- 1

while (sim_t < max_t) {
    print(sim_t)
    data_at_t <- metapop_data %>%
        filter(day == sim_t)
    
    if (sim_t <= 1) {
        data_at_t <- data_at_t %>%
            left_join(colonised_t, by = "fID") %>%
            mutate(susceptible = pmax(0, capacity_t - model_colonised_t))
    } else {
        data_at_t <- data_at_t %>%
            left_join(sim_results[[sim_t - dt]] %>% select(fID, model_colonised_t), by = "fID") %>%
            mutate(susceptible = pmax(0, capacity_t - model_colonised_t))
    }
    
    c_it <- if_else(
        data_at_t$capacity_t == 0,
        0,
        pmax(0, data_at_t$model_colonised_t +
            beta * data_at_t$susceptible * data_at_t$model_colonised_t / data_at_t$capacity_t +
            gamma * data_at_t$admissions_t -
            delta * data_at_t$model_colonised_t -
            (data_at_t$model_colonised_t / data_at_t$capacity_t) * data_at_t$discharges_t)
    )

    sim_results[[sim_t]] <- data_at_t %>% mutate(
        model_colonised_t = c_it
    )
    sim_t <- sim_t + dt
}

sim_results_bound <- bind_rows(sim_results)

max(sim_results_bound$model_colonised_t)

sim_results_bound %>%
    select(fID, day, model_colonised_t, capacity_t) %>%
    pivot_longer(!c(day, fID)) %>%
    ggplot(aes(x = day, y = value, colour = name)) +
    geom_line() +
    facet_wrap(~fID)

ggplot(sim_results_bound, aes(x=day, y=model_colonised_t, group=fID)) +
    geom_line()

ggplot(metapop_data, aes(x = day, y = capacity_t, group = fID)) +
    geom_line()
