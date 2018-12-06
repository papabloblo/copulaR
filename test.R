#' 
#' Carga de datos para las pruebas de desarrollo
#' 


# LIBRERÍAS ---------------------------------------------------------------

library(tidyverse)


# KDD98 -------------------------------------------------------------------

train <- read_csv("data/KDD1988_train.csv") %>% 
  mutate(id_data = "train")

test <- read_csv('data/KDD1988_test.csv') %>% 
  mutate(id_data = "test")


df <- bind_rows(train, test) %>% 
  # Elimar símbolo del dólar y convertir a numéricas
  mutate_at(c("TargetD",
              
              "GiftAvgLast",     
              "GiftAvg36",
              "GiftAvgAll", 
              
              "DemMedHomeValue",
              "DemMedIncome"
              ),
            
            function(x){ as.numeric(str_remove_all(x, "[\\$,]"))}
            )


train <- df %>% filter(id_data == "train") %>% select(-id_data, -ID)
target <- "TargetD"


valid = train
test = train
num_iter = 10
early_stopping_round = 0
num_sim = 500
max_bins = NULL
num_obs_fit = NULL
bin_target = FALSE
eval_metric  = "MAPE"
verbosity = TRUE