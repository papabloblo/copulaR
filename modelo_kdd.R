library(tidyverse)


# Carga de datos ----------------------------------------------------------
train_kdd <- read_csv("data/KDD1988_train.csv")
test_kdd <- read_csv("data/KDD1988_test.csv")

# Tratamiento de datos ----------------------------------------------------
train_kdd$id_tabla <- "train"
test_kdd$id_tabla <- "test"

kdd <- bind_rows(train_kdd, test_kdd)

kdd <- kdd %>%
  mutate(TargetD         = drop_dollar_comma(TargetD),
         GiftAvgLast     = drop_dollar_comma(GiftAvgLast),
         GiftAvg36       = drop_dollar_comma(GiftAvg36),
         GiftAvgAll      = drop_dollar_comma(GiftAvgAll),
         DemMedHomeValue = drop_dollar_comma(DemMedHomeValue),
         DemMedIncome    = drop_dollar_comma(DemMedIncome))

train_kdd <- kdd %>%
  filter(id_tabla == "train") %>%
  select(-id_tabla)

test_kdd <- kdd %>%
  filter(id_tabla == "test") %>%
  select(-id_tabla)
