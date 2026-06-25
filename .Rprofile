# Project startup code (run once at R startup, after .Renviron).
# -----------------------------------------------------------------------------
# Cap arrow's compute/IO threads so heavy work (e.g. the HES Parquet conversion)
# never saturates the shared server. This is a function call, not an env var, so
# it lives here rather than in .Renviron.
#
# It is set via arrow's onLoad hook rather than called directly: that way it only
# runs if and when arrow is actually loaded, so the session still starts cleanly
# in this project even if a given task does not use arrow. Adjust the numbers if
# the server's core budget changes (keep them in step with .Renviron).

setHook(
  packageEvent("arrow", "onLoad"),
  function(...) {
    try(arrow::set_cpu_count(6),       silent = TRUE)  # compute threads
    try(arrow::set_io_thread_count(4), silent = TRUE)  # IO / write threads
  }
)

# brief, non-intrusive note so an analyst knows the caps are active
if (interactive())
  message("Project .Rprofile: thread caps active (arrow 6 cpu / 4 io; OMP/BLAS 6).")
