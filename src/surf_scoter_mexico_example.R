# download only the season proportion of population layer
ebirdst_download_status("Surf Scoter", 
                        pattern = "proportion-population_seasonal_mean_3km")

# breeding season proportion of population
abd_nonbreeding <- load_raster("Surf Scoter",
                               product = "proportion-population",
                               period = "seasonal") |> 
  subset("nonbreeding")

# load a polygon for the boundary of Mexico
mexico <- ne_countries(country = "Mexico") |> 
  st_transform(crs = st_crs(abd_nonbreeding))

# proportion in mexico, wrong estimation because of coastal species
wrong_proportion <- extract(abd_nonbreeding, mexico, fun = "sum", na.rm = TRUE)
print(wrong_proportion)

### Proper estimated proportion buffering the coastal shore of Mexico

# buffer by 5000m = 5km
mexico_buffer <- st_buffer(mexico, dist = 5000)

# proportion in mexico
correct_proportion <- extract(abd_nonbreeding, mexico_buffer, fun = "sum", na.rm = TRUE,
        touches = TRUE)

print(correct_proportion
)
