# -----------------------------------------------------------------------------
# Load data (example)
# -----------------------------------------------------------------------------
dat = readRDS('.../data/...')

names(dat)
unique(dat$studyid)

studies <- c("WASH-Bangladesh", "WASH-Kenya")

setnames(dat, old = c("arm.x"), new = c ("arm"))

dat0 = filter(dat, studyid %in% studies) %>% 
  group_by(studyid, subjid) %>% 
  mutate_at("subjid",as.numeric) %>%
  arrange(studyid,subjid, agedays) %>%
  filter(row_number()==n()) 

# WASH-Kenya
wash_kenya <- dat0[dat0$studyid == "WASH-Kenya",]
wash_kenya <- wash_kenya %>%
  droplevels() %>%
  ungroup() %>%
  dplyr::select(agedays, haz, sex, month, arm, brthmon, cleanck, impfloor, enstunt,
                enwast,  W_mage, W_mhtcm, W_meducyrs, W_nhh, hhwealth_quart, impsan) %>%
  select_if(~!all(is.na(.))) %>%
  drop_na(arm, haz) %>%
  mutate(a = as.integer((arm != 'Control') & (arm != 'Passive Control'))) %>% select(-arm) %>% mutate(sex = as.factor(sex))

#summary(wash_kenya)

nodes0 <- list(W = colnames(wash_kenya) %w/o% c('haz', 'a'),
               A = 'a',
               Y = 'haz')

wash_kenya = process_missing(wash_kenya, nodes0)[[1]]

colnames(wash_kenya) <- make.names(colnames(wash_kenya), unique=TRUE)



# WASH-Bangladesh
wash_b <- dat0[dat0$studyid == "WASH-Bangladesh",]
wash_b <- wash_b %>% 
  droplevels() %>%
  ungroup() %>%
  dplyr::select(arm, haz,  W_mage, W_mhtcm, W_mwtkg, W_mbmi, sex, month,
                brthmon, hfoodsec, enstunt, enwast, impsan, hhwealth_quart, 
                agedays, W_meducyrs, W_feducyrs, W_nhh, W_parity) %>%
  select_if(~!all(is.na(.))) %>%
  drop_na(arm, haz) %>%
  mutate(a = as.integer(arm != 'Control')) %>% select(-arm) %>% mutate(sex = as.factor(sex))

#process missingness
nodes0 <- list(W = colnames(wash_b) %w/o% c('haz', 'a'),
               A = 'a',
               Y = 'haz')

wash_b = process_missing(wash_b, nodes0)[[1]]

colnames(wash_b) <- make.names(colnames(wash_b), unique=TRUE)


# save example data
save_path <- ".../data"
write.csv(wash_b, paste0(save_path, "/wash_b.csv"), row.names = FALSE)
write.csv(wash_kenya, paste0(save_path, "/wash_k.csv"), row.names = FALSE)
