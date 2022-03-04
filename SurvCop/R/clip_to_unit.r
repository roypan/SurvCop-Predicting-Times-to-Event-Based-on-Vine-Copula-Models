clip_to_unit <- function(value) {
  # clip value to (0, 1).
  pmax(pmin(value, 1 - 1e-7), 1e-7)
}