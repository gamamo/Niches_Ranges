# Load packages -----------------------------------------------------------
library(tidyverse)
library(arrow)
library(tictoc)
library(here)
library(fs)
library(terra)
library(tidyterra)


# Load the recent created dataset ------------------------------------------------------------
ds <- open_dataset(
  here("products","PAMs_global_moll_sum","ranges_ds"),
  format = "parquet")

dim(ds)

dir_tree(here("products","PAMs_global_moll_sum","ranges_ds"))

total_size <- sum(dir_info(
  here("products","PAMs_global_moll_sum","ranges_ds"), 
  recurse = TRUE)$size, na.rm = TRUE); total_size

ds |> head() |> collect() -> hd

