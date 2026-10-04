# =============================================================================
#  Section 3.9 - Is the second yellow card rate specific to La Liga?
#
#  This script replicates ONE quantity from Section 3.7 across all five major
#  European leagues: the probability that a booked player receives a second
#  yellow card. It then tests whether the five rates are distinguishable from
#  one another.
#
#  Scope, deliberately narrow. No expected goals model, no action value model,
#  no panel, no fixed effects, no coordinates. This matters: the coordinate
#  encoding differs between leagues (Italy, France and England store positions
#  as [{'y': 52, 'x': 49}] while Germany stores [(50,50)]), and an earlier
#  attempt to parse them uniformly silently produced 100% missing values in
#  three leagues. Card tags carry no such ambiguity, so this replication does
#  not touch the fragile part of the data.
#
#  Because we compare raw proportions rather than model estimates, there is
#  nothing to adjust for. League differences in refereeing, tempo or scoring
#  environment would matter if we were comparing conditional effects; they do
#  not enter a descriptive rate.
#
#  Run Analysis_secondYellow.R first. This script is standalone but the
#  Spain figures it reproduces are the validation gate for everything below.
# =============================================================================

library(tidyverse)

DATA_DIR  <- "data" # folder containing events_*.csv for all five leagues
MIN_AFTER <- 15      # minimum minutes remaining, as in Section 3.7
N_PERM    <- 10000   # permutation replicates
set.seed(20260914)   # so the permutation p-value is reproducible

# Wyscout card tags. These identifiers are shared across all leagues in the
# Pappalardo et al. (2019) release, unlike the positional encoding.
#
# Note on 1701 and 1702. An earlier version of this script had these two
# reversed. The error was caught by the parser check in extract_cards(), which
# compares the parsed flags against the independently produced columns in
# events_Spain.csv: the parsed yellow count came back as 29 and the parsed red
# count as 1863, exactly the reverse of the stored values. The assignment below
# is the one that reproduces the stored columns.
TAG_RED    <- 1701   # straight red
TAG_YELLOW <- 1702   # caution
TAG_SECOND <- 1703   # second caution, dismissal


# =============================================================================
# 0. Locate the event files
# =============================================================================
files <- list.files(DATA_DIR, pattern = "^events_.*\\.csv$", full.names = TRUE)
cat("========== 0. Files found ==========\n")
if (length(files) == 0) stop("No events_*.csv files found in DATA_DIR.")
for (f in files) cat(" ", basename(f), "\n")
cat("\n")


# =============================================================================
# 1. Helper functions
#
#    The two schemas in the data release are handled by one code path.
#    events_Spain.csv has already been processed and carries is_yellow,
#    is_second_yellow, is_red and match_time as columns. The other four are
#    raw and carry a tags string and eventSec/matchPeriod instead.
#
#    We deliberately parse the tags even for Spain, where the answer is known,
#    so that the parser can be checked against columns that were produced
#    independently. If the parser reproduces Spain exactly, the same parser
#    applied to the other four leagues can be trusted.
# =============================================================================

# Match a tag id that is not part of a longer number.
has_tag <- function(x, id) {
  grepl(paste0("(?<![0-9])", id, "(?![0-9])"), x, perl = TRUE)
}

# Absolute seconds from kickoff. Extra-time periods do not occur in league
# play but are mapped anyway so the function cannot silently return NA.
period_offset <- function(p) {
  dplyr::case_when(
    p == "1H" ~ 0,
    p == "2H" ~ 2700,
    p == "E1" ~ 5400,
    p == "E2" ~ 6300,
    p == "P"  ~ 7200,
    TRUE      ~ NA_real_
  )
}

extract_cards <- function(path) {
  ev <- read.csv(path, stringsAsFactors = FALSE)
  league <- sub("^events_", "", sub("\\.csv$", "", basename(path)))
  league <- sub("-.*$", "", league)   # strip any download suffix
  
  tag_col <- intersect(c("tags_raw", "tags"), names(ev))
  if (length(tag_col) == 0)
    stop("No tags column in ", basename(path),
         ". Columns present: ", paste(names(ev), collapse = ", "))
  tags <- ev[[tag_col[1]]]
  
  # Card flags, always from the tags, for every league.
  yellow <- has_tag(tags, TAG_YELLOW)
  second <- has_tag(tags, TAG_SECOND)
  red    <- has_tag(tags, TAG_RED)
  
  # Time, from the existing column where available and reconstructed otherwise.
  if ("match_time" %in% names(ev)) {
    mt <- ev$match_time
    time_source <- "match_time column"
  } else {
    mt <- ev$eventSec + period_offset(ev$matchPeriod)
    time_source <- "eventSec + period offset"
  }
  
  # ---- Parser validation, where independent columns exist to check against --
  check <- NULL
  if (all(c("is_yellow", "is_second_yellow", "is_red") %in% names(ev))) {
    check <- c(
      yellow = sum(yellow) - sum(as.logical(ev$is_yellow),        na.rm = TRUE),
      second = sum(second) - sum(as.logical(ev$is_second_yellow), na.rm = TRUE),
      red    = sum(red)    - sum(as.logical(ev$is_red),           na.rm = TRUE)
    )
  }
  
  # Any unrecognised matchPeriod would silently produce NA times, which would
  # then propagate into card_time and mins_left. Count them here instead.
  n_bad_time <- sum(is.na(mt))
  
  cat(sprintf("%-10s rows %7d | tags: %-9s | time: %-26s",
              league, nrow(ev), tag_col[1], time_source))
  cat(sprintf("| Y %5d  SY %3d  R %3d", sum(yellow), sum(second), sum(red)))
  if (n_bad_time > 0) cat(sprintf("  [%d rows with missing time]", n_bad_time))
  if (!is.null(check)) {
    if (all(check == 0)) cat("  [parser matches stored columns]\n")
    else {
      cat("  [MISMATCH vs stored columns: ",
          paste(names(check), check, collapse = ", "), "]\n")
      stop("Tag parser disagrees with the stored columns in ", basename(path),
           ". Fix the tag identifiers before going further.")
    }
  } else cat("\n")
  
  out <- data.frame(
    league     = league,
    playerId   = ev$playerId,
    matchId    = ev$matchId,
    match_time = mt,
    yellow     = yellow,
    second     = second,
    red        = red
  )
  rm(ev, tags, yellow, second, red, mt); gc(verbose = FALSE)
  out
}


# =============================================================================
# 2. Build the booked player-match table, one league at a time
#
#    Files are read and discarded individually rather than stacked, because
#    five full event streams held in memory at once is several gigabytes.
#
#    The logic below is copied from Section 2 of Analysis_secondYellow.R and
#    must stay identical to it, or the Spain row will not reproduce.
# =============================================================================
cat("========== 1. Reading and parsing ==========\n")

build_booked <- function(cards) {
  match_end <- cards |>
    group_by(matchId) |>
    summarise(match_end = max(match_time, na.rm = TRUE), .groups = "drop")
  
  first_yellow <- cards |>
    filter(yellow) |>
    group_by(playerId, matchId) |>
    summarise(card_time = min(match_time), .groups = "drop")
  
  dismissal <- cards |>
    filter(second) |>
    group_by(playerId, matchId) |>
    summarise(sy_time = min(match_time), .groups = "drop")
  
  first_yellow |>
    mutate(
      key        = paste(playerId, matchId),
      sy_time    = dismissal$sy_time[match(key, paste(dismissal$playerId, dismissal$matchId))],
      got_second = ifelse(is.na(sy_time), 0, 1),
      match_end  = match_end$match_end[match(matchId, match_end$matchId)]
    ) |>
    mutate(
      end_time  = ifelse(got_second == 1, sy_time, match_end),
      mins_left = (match_end - card_time) / 60
    ) |>
    filter(end_time > card_time) |>
    select(playerId, matchId, card_time, got_second, mins_left)
}

booked_all <- list()
n_matches  <- c()
for (f in files) {
  cards <- extract_cards(f)
  lg <- cards$league[1]
  n_matches[lg]   <- length(unique(cards$matchId))
  booked_all[[lg]] <- build_booked(cards) |> mutate(league = lg)
  rm(cards); gc(verbose = FALSE)
}
booked <- bind_rows(booked_all)
cat("\n")


# =============================================================================
# 3. Validation gate
#
#    Spain must reproduce Section 3.7 exactly: 1,365 booked player-matches with
#    at least fifteen minutes remaining, of which 38 ended in a dismissal, a
#    rate of 2.784%. If this row does not match, the pipeline has diverged from
#    the main analysis and no other league should be reported.
# =============================================================================
cat("========== 2. Validation gate (Spain) ==========\n")
sp <- booked |> filter(league == "Spain", mins_left >= MIN_AFTER)
cat("Matches         :", n_matches["Spain"], " (expected 380)\n")
cat("Booked, >=15 min:", nrow(sp),          " (expected 1365)\n")
cat("Dismissals      :", sum(sp$got_second)," (expected 38)\n")
cat("Rate            :", round(100*mean(sp$got_second), 3), "% (expected 2.784)\n")
gate_ok <- nrow(sp) == 1365 && sum(sp$got_second) == 38
if (!gate_ok) {
  cat("\nGATE FAILED.\n")
  stop("Spain does not reproduce Section 3.7. The pipeline has diverged from ",
       "the main analysis, so the other four leagues are not interpretable. ",
       "Stopping here rather than printing numbers that cannot be trusted.")
}
cat("\nGATE PASSED.\n\n")


# =============================================================================
# 4. Result - second yellow rate by league
# =============================================================================
cat("========== 3. Second yellow rate by league ==========\n")

analysis <- booked |> filter(mins_left >= MIN_AFTER)

by_league <- analysis |>
  group_by(league) |>
  summarise(booked = n(), dismissals = sum(got_second), .groups = "drop") |>
  rowwise() |>
  mutate(
    rate_pct = 100 * dismissals / booked,
    ci_lo    = 100 * binom.test(dismissals, booked)$conf.int[1],
    ci_hi    = 100 * binom.test(dismissals, booked)$conf.int[2]
  ) |>
  ungroup() |>
  arrange(rate_pct)

by_league |>
  mutate(across(c(rate_pct, ci_lo, ci_hi), ~ round(.x, 2))) |>
  as.data.frame() |> print()

pooled_n <- nrow(analysis); pooled_d <- sum(analysis$got_second)
pt <- binom.test(pooled_d, pooled_n)
cat(sprintf("\nPooled: %d bookings, %d dismissals, %.3f%%  95%% CI [%.2f, %.2f]\n",
            pooled_n, pooled_d, 100*pooled_d/pooled_n,
            100*pt$conf.int[1], 100*pt$conf.int[2]))

# Where does La Liga sit relative to the others?
sp_rate <- by_league$rate_pct[by_league$league == "Spain"]
cat(sprintf("La Liga rate %.3f%% ranks %d of %d leagues.\n",
            sp_rate, which(by_league$league == "Spain"), nrow(by_league)))


# =============================================================================
# 5. Are the five rates distinguishable?
#
#    Three tests of the same null hypothesis, that all five leagues share one
#    underlying rate. They use different machinery, so agreement between them
#    is informative and disagreement would be a warning.
#
#      (a) Pearson chi-square, relying on a large-sample approximation
#      (b) Monte Carlo Fisher exact test
#      (c) Permutation test, shuffling the dismissal indicator across leagues
#
#    The permutation test is included because the dismissal counts per league
#    are small (tens, not hundreds), which is exactly the regime where the
#    chi-square approximation is worth checking rather than assuming.
# =============================================================================
cat("\n========== 4. Test of homogeneity ==========\n")

# Built from the per-league summary rather than by reshaping, so that the two
# columns exist even if some league has no dismissals at all.
tab <- as.matrix(data.frame(
  survived  = by_league$booked - by_league$dismissals,
  dismissed = by_league$dismissals,
  row.names = by_league$league
))
cat("\nContingency table:\n"); print(tab)

cat("\nExpected dismissals under a common rate:\n")
p0 <- pooled_d / pooled_n
exp_tbl <- data.frame(
  league   = rownames(tab),
  observed = tab[, "dismissed"],
  expected = round(rowSums(tab) * p0, 1)
)
print(exp_tbl, row.names = FALSE)

# ---- (a) chi-square ---------------------------------------------------------
cs <- suppressWarnings(chisq.test(tab))
cat(sprintf("\n(a) Pearson chi-square: X2 = %.3f, df = %d, p = %.4f\n",
            cs$statistic, cs$parameter, cs$p.value))
cat("    Minimum expected cell count:", round(min(cs$expected), 1),
    ifelse(min(cs$expected) >= 5, "(approximation acceptable)",
           "(approximation questionable)"), "\n")

# ---- (b) Fisher, Monte Carlo ------------------------------------------------
fi <- fisher.test(tab, simulate.p.value = TRUE, B = N_PERM)
cat(sprintf("(b) Fisher exact (Monte Carlo, B = %d): p = %.4f\n", N_PERM, fi$p.value))

# ---- (c) permutation --------------------------------------------------------
#     League sizes and the total number of dismissals are both held fixed;
#     only the assignment of dismissals to leagues is randomized.
gidx  <- as.integer(factor(analysis$league))
n_l   <- as.vector(table(analysis$league))
exp_d <- n_l * p0
chi_fast <- function(d) {
  sum((d - exp_d)^2 / exp_d) + sum((d - exp_d)^2 / (n_l - exp_d))
}
obs_stat  <- chi_fast(as.vector(rowsum(analysis$got_second, gidx)))
perm_stat <- replicate(N_PERM, chi_fast(as.vector(rowsum(sample(analysis$got_second), gidx))))
p_perm    <- (1 + sum(perm_stat >= obs_stat)) / (1 + N_PERM)
cat(sprintf("(c) Permutation (B = %d): observed X2 = %.3f, p = %.4f\n",
            N_PERM, obs_stat, p_perm))

cat("\nInterpretation: a large p-value means the five leagues are consistent\n")
cat("with a single underlying rate, i.e. the La Liga figure used in Section\n")
cat("3.8 is not a property of that league in particular.\n")


# =============================================================================
# 6. Pairwise comparison - exploratory only
#
#    The widest contrast will look dramatic. It is reported here with a
#    Bonferroni correction so that the correction can be seen rather than
#    described, and it should not be quoted without it: ten comparisons drawn
#    from five leagues will produce an apparent extreme by construction.
# =============================================================================
cat("\n========== 5. Pairwise (exploratory) ==========\n")
lgs <- rownames(tab); pairs_out <- list()
for (i in 1:(length(lgs)-1)) for (j in (i+1):length(lgs)) {
  m <- tab[c(lgs[i], lgs[j]), c("dismissed","survived")]
  ft <- fisher.test(m)
  pairs_out[[length(pairs_out)+1]] <- data.frame(
    a = lgs[i], b = lgs[j], p_raw = ft$p.value)
}
pw <- bind_rows(pairs_out) |>
  mutate(p_bonferroni = pmin(1, p_raw * n())) |>
  arrange(p_raw)
pw |> mutate(across(starts_with("p_"), ~ round(.x, 4))) |> as.data.frame() |> print()
cat("\nComparisons made:", nrow(pw),
    " / surviving Bonferroni at 0.05:", sum(pw$p_bonferroni < 0.05), "\n")


# =============================================================================
# 7. Figure 9
# =============================================================================
fig <- by_league |>
  mutate(league = factor(league, levels = by_league$league)) |>
  ggplot(aes(x = league, y = rate_pct)) +
  geom_hline(yintercept = 100*p0, linetype = "dashed", colour = "grey40") +
  geom_errorbar(aes(ymin = ci_lo, ymax = ci_hi), width = 0.12, colour = "grey35") +
  geom_point(size = 3, colour = "firebrick") +
  labs(x = NULL, y = "Second yellow rate (%)",
       title = "Probability of a second yellow card, given a first",
       subtitle = paste0("2017-18 seasons; bookings with at least ", MIN_AFTER,
                         " minutes remaining. Dashed line is the pooled rate.")) +
  theme_minimal(base_size = 13)

print(fig)
ggsave(file.path(DATA_DIR, "fig_cross_league_sy.png"), fig,
       width = 8, height = 4.5, dpi = 300, bg = "white")

cat("\nDone. Figure written to fig_cross_league_sy.png in the data directory.\n")

cat("\n========== Session information ==========\n")
print(sessionInfo())