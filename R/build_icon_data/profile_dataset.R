# profile_dataset.R
#
# Describes a data frame you have not looked at column by column yet. For each
# column it works out what sort of thing it holds, then prints a summary that
# suits it: dates get a range, measurements get an average and quantiles,
# categories get a count per group, and anything that looks like a patient
# identifier or free text is only described by its shape, never its contents.
#
# Usage:
#   source("profile_dataset.R")
#   profile_dataset(cwt_og)
#
# To also save the output to a text file, alongside printing it:
#   profile_dataset(cwt_og, logfile = "cwt_og_profile.txt")
#
# This is the R equivalent of profile_cohort.do, for looking at a data frame
# already loaded in R rather than a .dta file on disk.

profile_dataset <- function(df, hidesmall = 10, roundto = 5, showextremes = FALSE,
                            maxgroups = 200, numgroups = 25, dateparse = 0.90,
                            iduniq = 0.98, logfile = NULL) {
  
  out <- character(0)                 # everything we print gets collected here too,
  # so it can be written to a log file at the end
  say <- function(...) {
    line <- paste0(...)
    cat(line, "\n", sep = "")
    out <<- c(out, line)
  }
  
  nm <- deparse(substitute(df))
  N <- nrow(df)
  K <- ncol(df)
  
  say(strrep("-", 78))
  say("Looking at: ", nm)
  say("Date: ", format(Sys.time(), "%d %b %Y %H:%M"))
  if (hidesmall > 1 || roundto > 1) {
    say("Small groups hidden, counts rounded to the nearest ", roundto, ".")
  } else {
    say("Exact counts. Identifiers and long notes are not printed.")
  }
  say("This never shows individual rows, patient ids, or free text.")
  say("It only ever shows one variable at a time, never a cross-tab of two")
  say("variables together, which keeps the risk of identifying anyone much")
  say("lower than a table would carry.")
  say("Even so, treat this like any other output: run it past your usual")
  say("checks before anything from it leaves your organisation.")
  say(strrep("-", 78))
  say("")
  say("Rows: ", N, "   Columns: ", K)
  
  if (N == 0) {
    say("This data frame has no rows in it.")
    return(invisible(NULL))
  }
  
  # name-based hints, checked once per column
  id_patterns   <- c("patientid", "pseudo", "tumourid", "nhsnumber", "lsoa", "postcode", "_id")
  org_patterns  <- c("trust", "hosp", "provider", "site", "org", "procode", "alliance", "ccg")
  code_patterns <- c("morph", "icd_", "icdcode", "_icd", "opcs", "topog")
  
  name_hints <- function(varname) {
    ln <- tolower(varname)
    nameid  <- any(vapply(id_patterns, function(p) grepl(p, ln, fixed = TRUE), logical(1))) || ln == "id"
    nameorg <- any(vapply(org_patterns, function(p) grepl(p, ln, fixed = TRUE), logical(1)))
    namecode <- any(vapply(code_patterns, function(p) grepl(p, ln, fixed = TRUE), logical(1))) || ln == "type"
    list(nameid = nameid, nameorg = nameorg | namecode)
  }
  
  # tries a handful of date orders on a text column and keeps whichever reads
  # the most values as a real date. returns NULL if none of them read well.
  try_parse_dates <- function(x) {
    xx <- x[!is.na(x) & x != ""]
    if (length(xx) == 0) return(NULL)
    orders <- c("%d/%m/%Y", "%m/%d/%Y", "%Y-%m-%d", "%d-%m-%Y", "%Y/%m/%d")
    best <- NULL
    best_ok <- 0
    for (fmt in orders) {
      parsed <- suppressWarnings(as.Date(xx, format = fmt))
      ok <- sum(!is.na(parsed))
      if (ok > best_ok) {
        best_ok <- ok
        best <- list(fmt = fmt, ok = ok)
      }
    }
    if (is.null(best) || best$ok / length(xx) < dateparse) return(NULL)
    best
  }
  
  # rounds a set of counts and works out which, if any, should be hidden for
  # being too small a group. if exactly one group would be hidden, a second
  # (the next smallest) is hidden too, so it cannot be worked out by taking
  # the others away from the total.
  suppress_counts <- function(counts) {
    hide <- rep(FALSE, length(counts))
    if (hidesmall > 1) {
      hide <- counts < hidesmall
      if (sum(hide) == 1 && length(counts) > 1) {
        visible <- which(!hide)
        smallest <- visible[which.min(counts[visible])]
        hide[smallest] <- TRUE
      }
    }
    hide
  }
  
  print_groups <- function(labels, counts, nfilled, always_list = FALSE) {
    ngroups <- length(labels)
    if (ngroups > maxgroups && !always_list) {
      say("    ", ngroups, " different groups, too many to list here.")
      return(invisible(NULL))
    }
    hide <- suppress_counts(counts)
    base <- if (roundto > 1) round(nfilled / roundto) * roundto else nfilled
    
    if (hidesmall > 1 || roundto > 1) {
      say("    groups (out of ", base, " rows with a value; small groups hidden, counts rounded):")
      if (sum(hide) > sum(counts < hidesmall)) {
        say("      [a second group is hidden so the first cannot be worked out by subtraction]")
      }
    } else {
      say("    groups (out of ", base, " rows with a value):")
    }
    
    ord <- order(labels)
    for (i in ord) {
      if (hide[i]) {
        say("      ", labels[i], "   n=hidden   (.%)")
      } else {
        shown_n <- if (roundto > 1) round(counts[i] / roundto) * roundto else counts[i]
        pct <- if (base > 0) sprintf("%.1f", 100 * shown_n / base) else "0.0"
        say("      ", labels[i], "   n=", shown_n, "   (", pct, "%)")
      }
    }
  }
  
  summarise_measurement <- function(x) {
    xx <- x[!is.na(x)]
    qs <- quantile(xx, probs = c(.01, .05, .25, .5, .75, .95, .99), names = FALSE, type = 7)
    say(sprintf("    average %.2f   spread %.2f", mean(xx), sd(xx)))
    if (showextremes) {
      say("    smallest ", format(min(xx)), "   largest ", format(max(xx)))
    }
    say("    1% / 5% / 25% / half / 75% / 95% / 99% of rows fall below:")
    say("      ", paste(signif(qs, 5), collapse = " / "))
    lopsided <- mean(xx) - qs[4]
    say(sprintf("    (average minus middle value = %.2f. a positive figure means a few large values pull the average up)", lopsided))
    nneg <- sum(xx < 0)
    if (nneg > 0) say("    note: ", nneg, " rows hold a negative value, worth a look")
  }
  
  # ---- go through every column -------------------------------------------
  say("")
  say(strrep("-", 78))
  say("1. Each column in turn")
  say(strrep("-", 78))
  
  manifest <- data.frame(column = character(0), class = character(0),
                         ndiff = integer(0), pblank = numeric(0),
                         lookup = character(0), stringsAsFactors = FALSE)
  
  for (v in names(df)) {
    x <- df[[v]]
    hints <- name_hints(v)
    nameid <- hints$nameid
    nameorg <- hints$nameorg
    
    nblank <- sum(is.na(x)) + if (is.character(x)) sum(x == "" & !is.na(x)) else 0
    nfilled <- N - nblank
    pblank <- 100 * nblank / N
    
    lookup_used <- ""
    class <- NA_character_
    ndiff <- NA_integer_
    
    say("")
    
    # -- dates: proper R Date/POSIXct columns ------------------------------
    if (inherits(x, "Date") || inherits(x, "POSIXct")) {
      class <- "date"
      xx <- as.Date(x[!is.na(x)])
      ndiff <- length(unique(xx))
      say(v, "  [", class(x)[1], ", date]")
      say("    blank: ", nblank, " (", sprintf("%.1f", pblank), "%)   different values: ", ndiff)
      if (length(xx) > 0) {
        say("    runs from ", format(min(xx), "%d %b %Y"), " to ", format(max(xx), "%d %b %Y"),
            ", middle date ", format(median(xx), "%d %b %Y"))
      } else {
        say("    (every row is blank, so nothing to summarise)")
      }
      
      # -- factor columns: already a set of categories -----------------------
    } else if (is.factor(x)) {
      class <- "groups"
      lv <- levels(x)
      xna <- x[!is.na(x)]
      ndiff <- length(unique(xna))
      say(v, "  [factor, groups]")
      say("    blank: ", nblank, " (", sprintf("%.1f", pblank), "%)   different values: ", ndiff)
      if (nfilled == 0) {
        say("    (every row is blank, so nothing to summarise)")
      } else {
        tab <- table(xna)
        print_groups(names(tab), as.integer(tab), nfilled, always_list = nameorg)
      }
      
      # -- labelled numerics (haven-style, carrying a value label) -----------
    } else if (inherits(x, "haven_labelled")) {
      class <- "groups"
      labs <- attr(x, "labels")
      lookup_used <- paste(names(labs), collapse = ",")
      xna <- x[!is.na(x)]
      ndiff <- length(unique(xna))
      say(v, "  [labelled, groups]")
      vlab <- attr(x, "label", exact = TRUE)
      if (!is.null(vlab) && nzchar(vlab)) say("    description: ", vlab)
      say("    blank: ", nblank, " (", sprintf("%.1f", pblank), "%)   different values: ", ndiff)
      if (nfilled == 0) {
        say("    (every row is blank, so nothing to summarise)")
      } else {
        tab <- table(xna)
        codes <- as.numeric(names(tab))
        shown <- vapply(codes, function(cd) {
          hit <- which(labs == cd)
          if (length(hit) == 0) paste0(cd, " = (no meaning attached to this code)")
          else paste0(cd, " = ", names(labs)[hit[1]])
        }, character(1))
        print_groups(shown, as.integer(tab), nfilled, always_list = nameorg)
        # unused levels: defined in the label but never occurring in the data
        unused <- labs[!(labs %in% codes)]
        if (length(unused) > 0) {
          say("    (this column's lookup also defines, but never uses: ",
              paste(sprintf("%s=%s", unused, names(unused)), collapse = ", "), ")")
        }
      }
      
      # -- plain character columns --------------------------------------------
    } else if (is.character(x)) {
      xna <- x[!is.na(x) & x != ""]
      ndiff <- length(unique(xna))
      uniqshare <- if (N > 0) ndiff / N else 0
      manyvalues <- ndiff > numgroups
      typicallen <- if (length(xna) > 0) median(nchar(xna)) else 0
      looksnotes <- (!nameorg) && ((typicallen > 40 && ndiff > 50) || (uniqshare > 0.25 && ndiff > 100))
      dateguess <- if (nfilled > 0) try_parse_dates(xna) else NULL
      
      if (!is.null(dateguess)) {
        class <- "date"
        say(v, "  [character, date]")
        say("    blank: ", nblank, " (", sprintf("%.1f", pblank), "%)   different values: ", ndiff)
        parsed <- suppressWarnings(as.Date(xna, format = dateguess$fmt))
        parsed <- parsed[!is.na(parsed)]
        say("    stored as text, read using the pattern ", dateguess$fmt, ".")
        say("    runs from ", format(min(parsed), "%d %b %Y"), " to ", format(max(parsed), "%d %b %Y"),
            ", middle date ", format(median(parsed), "%d %b %Y"))
        bad <- dateguess$ok < length(xna)
        if (bad) say("    note: ", length(xna) - dateguess$ok, " rows could not be read as a date")
      } else if (nameid && manyvalues) {
        class <- "identifier"
        say(v, "  [character, identifier]")
        say("    blank: ", nblank, " (", sprintf("%.1f", pblank), "%)   different values: ", ndiff)
        say("    contents not printed. typical length ", round(typicallen), " characters")
      } else if (!nameorg && manyvalues && uniqshare > 0.5) {
        class <- "identifier"
        say(v, "  [character, identifier]")
        say("    blank: ", nblank, " (", sprintf("%.1f", pblank), "%)   different values: ", ndiff)
        say("    contents not printed. typical length ", round(typicallen), " characters")
      } else if (looksnotes) {
        class <- "notes"
        say(v, "  [character, notes]")
        say("    blank: ", nblank, " (", sprintf("%.1f", pblank), "%)   different values: ", ndiff)
        say("    contents not printed. typical length ", round(typicallen), " characters")
      } else {
        class <- "groups"
        say(v, "  [character, groups]")
        say("    blank: ", nblank, " (", sprintf("%.1f", pblank), "%)   different values: ", ndiff)
        if (nfilled == 0) {
          say("    (every row is blank, so nothing to summarise)")
        } else {
          tab <- table(xna)
          print_groups(names(tab), as.integer(tab), nfilled, always_list = nameorg)
        }
      }
      
      # -- plain numeric columns ----------------------------------------------
    } else if (is.numeric(x)) {
      xna <- x[!is.na(x)]
      ndiff <- length(unique(xna))
      uniqshare <- if (N > 0) ndiff / N else 0
      manyvalues <- ndiff > numgroups
      wholenumbers <- length(xna) > 0 && all(abs(xna - round(xna)) < 1e-9)
      ln <- tolower(v)
      namedate <- grepl("date", ln) || grepl("mdy", ln)
      
      if (namedate) {
        class <- "date"
        say(v, "  [", class(x)[1], ", date]")
        say("    blank: ", nblank, " (", sprintf("%.1f", pblank), "%)   different values: ", ndiff)
        if (nfilled == 0) {
          say("    (every row is blank, so nothing to summarise)")
        } else {
          say("    named like a date but held as a plain number. middle value ", median(xna))
        }
      } else if (nameid && (manyvalues || uniqshare >= 0.99)) {
        class <- "identifier"
        say(v, "  [", class(x)[1], ", identifier]")
        say("    blank: ", nblank, " (", sprintf("%.1f", pblank), "%)   different values: ", ndiff)
        say("    a number used as an identifier. contents not printed.")
      } else if (!nameorg && wholenumbers && manyvalues && uniqshare > iduniq) {
        class <- "identifier"
        say(v, "  [", class(x)[1], ", identifier]")
        say("    blank: ", nblank, " (", sprintf("%.1f", pblank), "%)   different values: ", ndiff)
        say("    a number used as an identifier. contents not printed.")
      } else if (!nameorg && manyvalues) {
        class <- "measurement"
        say(v, "  [", class(x)[1], ", measurement]")
        say("    blank: ", nblank, " (", sprintf("%.1f", pblank), "%)   different values: ", ndiff)
        if (nfilled == 0) {
          say("    (every row is blank, so nothing to summarise)")
        } else {
          summarise_measurement(xna)
        }
      } else {
        class <- "groups"
        say(v, "  [", class(x)[1], ", groups]")
        say("    blank: ", nblank, " (", sprintf("%.1f", pblank), "%)   different values: ", ndiff)
        if (nfilled == 0) {
          say("    (every row is blank, so nothing to summarise)")
        } else {
          tab <- table(xna)
          print_groups(names(tab), as.integer(tab), nfilled, always_list = nameorg)
        }
      }
      
      # -- logical columns: a two-group TRUE/FALSE tabulation ----------------
    } else if (is.logical(x)) {
      class <- "groups"
      xna <- x[!is.na(x)]
      ndiff <- length(unique(xna))
      say(v, "  [logical, groups]")
      say("    blank: ", nblank, " (", sprintf("%.1f", pblank), "%)   different values: ", ndiff)
      if (nfilled == 0) {
        say("    (every row is blank, so nothing to summarise)")
      } else {
        tab <- table(xna)
        print_groups(names(tab), as.integer(tab), nfilled)
      }
      
      # -- anything else: report the type and move on -------------------------
    } else {
      class <- "unclassified"
      say(v, "  [", class(x)[1], ", not recognised - skipped]")
    }
    
    manifest <- rbind(manifest, data.frame(column = v, class = class, ndiff = ifelse(is.na(ndiff), 0, ndiff),
                                           pblank = round(pblank, 1), lookup = lookup_used,
                                           stringsAsFactors = FALSE))
  }
  
  # ---- summary table --------------------------------------------------------
  say("")
  say(strrep("-", 92))
  say("2. Summary table (column | what it is | different values | % blank)")
  say(strrep("-", 92))
  for (i in seq_len(nrow(manifest))) {
    say(sprintf("  %-30s | %-11s | %6d | %5.1f%%",
                manifest$column[i], manifest$class[i], manifest$ndiff[i], manifest$pblank[i]))
  }
  
  say("")
  say(strrep("-", 78))
  say("Finished. No individual rows were printed. Identifiers and long notes")
  say("were left out on purpose. Do glance down the summary table in case")
  say("anything has ended up in the wrong category.")
  say(strrep("-", 78))
  
  if (!is.null(logfile)) {
    writeLines(out, logfile)
    cat("\n(also written to ", logfile, ")\n", sep = "")
  }
  
  invisible(manifest)
}