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
  here("ranges_ds"),
  format = "parquet")

dim(ds)

dir_tree(here("ranges_ds"))

total_size <- sum(dir_info(
  here("ranges_ds"), 
  recurse = TRUE)$size, na.rm = TRUE); total_size

ds |> head() |> collect() -> hd


