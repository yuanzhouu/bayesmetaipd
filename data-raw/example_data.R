## Code to prepare `example_data` directly from official Simulation Study 1

devtools::load_all()

example_data <- load_example()

ipd <- example_data$ipd
ad_nested <- example_data$ad_nested
ad_subgroup <- example_data$ad_subgroup
ad_partial <- example_data$ad_partial
subgroup <- example_data$subgroup
formula <- example_data$formula
nested_formula <- example_data$nested_formula
drm_formula <- example_data$drm_formula
nested_reported <- example_data$nested_reported
partial_terms <- example_data$partial_terms

save(example_data, ipd, ad_nested, ad_subgroup, ad_partial, subgroup,
     formula, nested_formula, drm_formula, nested_reported, partial_terms,
     file = "data/example_data.rda", compress = "xz")
cat("example_data.rda regenerated successfully!\n")
