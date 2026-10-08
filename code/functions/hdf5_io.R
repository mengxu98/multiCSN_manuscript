read_hdf5_category <- function(path, field) {
  codes <- rhdf5::h5read(path, paste0("/obs/", field))
  levels <- rhdf5::h5read(path, paste0("/obs/__categories/", field))
  levels[codes + 1L]
}
