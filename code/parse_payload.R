###------------------------------------###
###---USER-MADE VARIABLE NAMES---------###
###------------------------------------###
# Edit the following variables to tailor
# code to specific redcap project

study_arm = "rand_cohort" # column name indicating study arm
randomization_ready = "in_rand___1"
# baseline_event = "Baseline"
bal_covariates = c(
  # List covariate column names used to
  # compute balance
  "rand_race", "rand_gi", "rand_prep"
)
# study_center = "center"

min_n_algorithm <- 5 # Sample size required before starting adaptive randomization

###------------------------------------###
###------------------------------------###


###---------------------------###
### Libraries
###---------------------------###
library(dplyr)
library(tidyr)
library(redcapAPI)
library(RCurl)


###---------------------------###
### Get arguments from payload
###---------------------------###
args <- commandArgs(trailingOnly = TRUE)
pl <- data.frame(pid = args[1], record = args[2], instrument = args[3])

# write.csv(pl, "args.csv")



###---------------------------###
### Read in REDCap data
###---------------------------###
# use Sys.setevn(TOKEN = *REDCAP API TOKEN*)
token = Sys.getenv("TOKEN")
URL <- "https://redcap.nubic.northwestern.edu/redcap/api/"

con <- redcapAPI::redcapConnection(
  url = URL,
  token = token
)
data <- redcapAPI::exportRecords(con)

saveRDS(data, paste0("data/redcap_db_", Sys.Date(), ".RDS"))




###---------------------------###
### Separate by new/old units
###---------------------------###
###---New units are those that trigger the DET
new_obs <- data %>%
  filter(record_id %in% pl$record) 

old_obs <- data %>%
  filter( # Get data that is neither the person nor their partner
    !(record_id %in% pl$record) &
      !(record_id_p %in% pl$record))

write.csv(new_obs, paste0("data/new_obs_", Sys.Date(), ".csv"))
# write.csv(old_obs, paste0("old_obs_", Sys.Date(), ".csv"))

###---------------------------###
### Check if randomization
### needs to occur
###---------------------------###
new_to_randomize <- new_obs %>%
  filter(
    is.na(get(study_arm)),
    get(randomization_ready) == "Checked"
  )

# Determine if we are randomizing a single participant or a couple
single <- new_to_randomize$enrollment_group == "Dyad (Enrolled as Couple)"

if(!single){
  partner_data <- data %>% 
    filter(record_id %in% new_to_randomize$record_id_p)
}

already_randomized <- data %>%
  filter(!is.na(get(study_arm)))

write.csv(new_to_randomize, "data/ntr.csv", row.names = F)
write.csv(already_randomized, "data/ar.csv", row.names = F)

if(nrow(new_to_randomize) > 0){
  write(
    paste(
      "There are",
      nrow(new_to_randomize),
      "units to randomize. Running randomization script."
    ),
    "out.txt"
  )
  # write(getwd(), "wd.txt")
  source("randomization.R")
  
} else {
  
  write(
    paste(
      "There are no eligible units to randomize:\n",
      new_obs %>% pull(get(study_arm)) %>% is.na() %>% sum(),
      "recently modified record has yet to be randomized.\n",
      new_obs %>% pull(get(randomization_ready)) %>% is.na() %>% sum(),
      "recently modified records are indicated as 'Ready to randomize'."
    ),
    "out.txt")
}
