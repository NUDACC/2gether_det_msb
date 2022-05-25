source("source_msb.R")

## Notes
# This script is triggered by parse_payload.R
# parse_payload.R passes the following arguments:
#   new_to_randomize: Tibble of baseline data on units to randomize
#   already_randomized: Tibble of baseline data on units already randomized
#   bal_covariates: list of strings indicating variable names for balance
#   study_center: string indicating variable name for center.
#   study_arm: string indicating variable name for study arm.
#   min_n_algorithm: interger Sample size required before starting adaptive randomization


###---------------------------###
### Randomize
###---------------------------###
if(nrow(already_randomized) == 0){already_randomized = tibble()}

msb_output <- get_votes(
  data = already_randomized,
  new_data = new_to_randomize,
  # center = study_center,
  single = single,
  covariates = bal_covariates,
  treatment = study_arm,
  min_n_adapt = min_n_algorithm, 
  show_votes = T
)

new_arm <- generate_assignment(msb_output$prob)
rand_note <- paste(
  "Record", pl$record, 
  "assignment = ", new_arm, 
  "probability =", msb_output$prob, 
  "votes =", unlist(msb_output$votes), 
  "majority =", msb_output$majority, 
  Sys.Date())
write(rand_note, "new_randomization.txt", append = TRUE)

###---------------------------###
### Update REDCap
###---------------------------###

if(single){
  ntr <- new_to_randomize %>% 
    mutate(
      !!study_arm := new_arm, 
      randomization_notes = rand_note,
      rand_prob = msb_output$prob, 
      rand_votes = unlist(msb_output$votes)
      )
} else {
  # Couples get randomized together
  
  ntr <- bind_rows(new_to_randomize, partner_data) %>% 
    mutate(
      !!study_arm := new_arm, 
      randomization_notes = rand_note,
      rand_prob = msb_output$prob, 
      rand_votes = unlist(msb_output$votes)
    )
  
}
write.csv(ntr %>% mutate(today = Sys.Date()), "data/randomized_patient.csv", append = TRUE)


rc_con <- redcapConnection(
  url = "https://redcap.nubic.northwestern.edu/redcap/api/",
  token = token, # Set API token as global variable
  conn = con,
  project = "7091",
  config = httr::config()
)

importRecords(
  rc_con,
  data = ntr,
  overwriteBehavior = "normal",
  returnContent = "count"
)

write(
  paste(
    "updated RCDB with a new randomization.",
    "\nRecord:", ntr$record_id
  ),
  "success.txt"
)
