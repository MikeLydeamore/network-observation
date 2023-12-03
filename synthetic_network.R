library(HospitalNetwork)

# Generate a synthetic network
# sID = Patient
# fID = Facility
# Adate = Admission date
# Dddate = Discharge date
n_subjects <- 10000
synthetic_db <- HospitalNetwork::create_fake_subjectDB_clustered(
    n_subjects = n_subjects,
    n_facilities = 100,
    avg_n_stays = 3,
    n_clusters = 10
)

# Patients who have an infection:

num_infected_patients <- round(n_subjects * 0.05)
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
write.csv(uniform_detection, "data/uniform_detection.csv", row.names=F)
