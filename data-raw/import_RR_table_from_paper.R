#RR table in paper was generated here:
#  C:\Users\jmend\Documents\JPSM\Dissertation\R packages\jmpaper3\examples-long\PEVS_ADM_2016
#See file:///C:/Users/jmend/Documents/JPSM/Dissertation/R%20packages/jmpaper3/examples-long/PEVS_ADM_2016/ADM_11_new_setup_20230424.html
#  Environment including RR table saved to env11_pevs_adm.rds

##C:/Users/jmend/Documents/JPSM/Dissertation/Writeups/Paper 3 - Nonresponse and design/Math/stratified domain calcs/domain_est_stratified_recodes.R"

library(magrittr)
library(dplyr)

env11_pevs_adm <- readRDS("C:/Users/jmend/Documents/JPSM/Dissertation/R packages/samplingNR/data-raw/env11_pevs_adm.rds")

pevs_adm_2016_rrs <-
  env11_pevs_adm$fvap_alloc50_table_for_output_with_total[1:91,] %>%
  select(h:rr_h) %>%
  rename(service = eps_svc,
         paygrade = eps_pg,
         age = eps_age,
         region = eps_region,
         sex = eps_sex,
         Nhat_h = N_h,
         rr_h = rr_h) %>%
  mutate(h = as.character(h))

#usethis::use_data(pevs_adm_2016_rr, overwrite = TRUE)
usethis::use_data(pevs_adm_2016_rrs)
