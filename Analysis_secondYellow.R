# =============================================================================
#  Second Yellow Card Risk and the Cost of Dismissal
#  La Liga 2017-18, Wyscout event data (Pappalardo et al. 2019)
#
#  This script produces, in order:
#    (1) P(second yellow | already booked)          -> Section 3.7 of the paper
#    (2) Cost of playing a man down, in team xG difference per 90
#    (3) The two combined into an expected cost and a substitution threshold
#    (4) Robustness checks and the starter-to-replacement value gap
#
#  Every number reported in Sections 3.7 to 3.8 of the manuscript comes from
#  this file. Run top to bottom; no other script is required.
#
#  Terminology note: we use "xG difference" throughout (xG for minus xG
#  against). This is a difference, not a differential.
# =============================================================================

library(tidyverse)
library(fixest)

# ---- Parameters -------------------------------------------------------------
DATA_DIR  <- "data" # set this to the folder containing events_Spain.csv and players.csv
MIN_AFTER <- 15    # minimum minutes remaining after a booking for the
# substitution decision to be considered live
MIN_BIN_N <- 20    # minimum observations for a descriptive point in Figure 8

events  <- read.csv(file.path(DATA_DIR, "events_Spain.csv"), stringsAsFactors = FALSE)
players <- read.csv(file.path(DATA_DIR, "players.csv"),      stringsAsFactors = FALSE)


# =============================================================================
# 0. Validation gates on the raw file
#
#    These four counts are cross-checked against the figures reported in the
#    original competition deck. If any of them moves, something upstream has
#    changed and nothing below should be trusted.
# =============================================================================
cat("\n========== 0. Raw data validation ==========\n")
cat("Rows:", nrow(events), "\n")
cat("Columns:", paste(names(events), collapse = ", "), "\n")

events <- events |>
  mutate(
    is_yellow        = as.logical(is_yellow),
    is_second_yellow = as.logical(is_second_yellow),
    is_red           = as.logical(is_red),
    duel_won         = as.logical(duel_won),
    is_goal          = as.logical(is_goal)
  )

cat("\nFirst yellows   :", sum(events$is_yellow,        na.rm = TRUE), " (expected 1863)\n")
cat("Second yellows  :", sum(events$is_second_yellow, na.rm = TRUE), " (expected 42)\n")
cat("Straight reds   :", sum(events$is_red,           na.rm = TRUE), " (expected 29)\n")
cat("Matches         :", length(unique(events$matchId)), " (expected 380)\n")


# =============================================================================
# 1. Expected goals model
#
#    Four predictors: distance to goal, distance squared, shot angle, and a
#    free-kick indicator. The angle is the angle subtended at the shot location
#    by the two goalposts, obtained from the law of cosines. It is included
#    because two shots at equal distance can face very different amounts of
#    open goal.
#
#    A linear probability model is used rather than logistic regression because
#    the quantity we need downstream is an additive expectation: these
#    probabilities are summed within 5-minute bins, and calibration on the
#    natural scale matters more than discrimination.
# =============================================================================
cat("\n========== 1. Expected goals model ==========\n")

PITCH_L <- 105; PITCH_W <- 68; GOAL_W <- 7.32
POST_L  <- PITCH_W/2 - GOAL_W/2
POST_R  <- PITCH_W/2 + GOAL_W/2

shot_angle <- function(x_m, y_m) {
  a <- sqrt((PITCH_L - x_m)^2 + (POST_L - y_m)^2)   # distance to near post
  b <- sqrt((PITCH_L - x_m)^2 + (POST_R - y_m)^2)   # distance to far post
  cosv <- (a^2 + b^2 - GOAL_W^2) / (2*a*b)          # law of cosines
  acos(pmin(pmax(cosv, -1), 1))                     # clamp for floating point
}

shots <- events |>
  filter(eventName == "Shot" |
           (eventName == "Free Kick" & subEventName == "Free kick shot")) |>
  filter(!is.na(start_x)) |>
  mutate(
    x_m       = start_x * PITCH_L/100,   # Wyscout coordinates are 0-100
    y_m       = start_y * PITCH_W/100,
    distance  = sqrt((PITCH_L - x_m)^2 + (PITCH_W/2 - y_m)^2),
    angle     = shot_angle(x_m, y_m),
    free_kick = as.numeric(eventName == "Free Kick"),
    goal      = as.numeric(is_goal)
  )

xg_model <- lm(goal ~ distance + I(distance^2) + angle + free_kick, data = shots)
print(summary(xg_model)$coefficients)

shots <- shots |>
  mutate(xG = pmin(pmax(predict(xg_model, newdata = shots), 0), 1))

# Calibration gate: summed predictions should recover the observed goal total.
cat("\nTotal predicted xG:", round(sum(shots$xG), 1),
    " / actual goals:", sum(shots$goal),
    " / ratio:", round(sum(shots$xG)/sum(shots$goal), 3), " (should be near 1.00)\n")


# =============================================================================
# 2. Timing of first bookings and dismissals
#
#    For each booked player-match we record when the caution was issued, when
#    (if ever) the second yellow followed, and when the match ended.
#
#    Two distinct durations are needed and they are easy to confuse:
#      mins_left  = time from the card to the end of the match. This is the
#                   window over which the substitution decision applies.
#      expo_mins  = time from the card until the player left the pitch, i.e.
#                   the window over which he was actually exposed to the risk
#                   of a second booking. For a player who is never dismissed
#                   the two are identical.
#
#    Distinguishing them matters. An earlier version of this analysis compared
#    post-card foul rates without doing so and found an apparent 3.3x increase
#    for dismissed players. That was an artefact: the second yellow is itself a
#    foul, and it terminates the observation window, so the rate is inflated by
#    construction. Under a fixed window the difference disappeared.
# =============================================================================
match_end <- events |>
  group_by(matchId) |>
  summarise(match_end = max(match_time), .groups = "drop")

first_yellow <- events |>
  filter(is_yellow) |>
  group_by(playerId, matchId) |>
  summarise(card_time = min(match_time), .groups = "drop")

dismissal <- events |>
  filter(is_second_yellow) |>
  group_by(playerId, matchId) |>
  summarise(sy_time = min(match_time), .groups = "drop")

booked <- first_yellow |>
  mutate(
    key        = paste(playerId, matchId),
    sy_time    = dismissal$sy_time[match(key, paste(dismissal$playerId, dismissal$matchId))],
    got_second = ifelse(is.na(sy_time), 0, 1),
    match_end  = match_end$match_end[match(matchId, match_end$matchId)]
  ) |>
  mutate(
    end_time  = ifelse(got_second == 1, sy_time, match_end),
    mins_left = (match_end - card_time) / 60,
    expo_mins = (end_time  - card_time) / 60
  ) |>
  filter(end_time > card_time)   # drops cards issued on the final recorded event


# =============================================================================
# 3. Result 1 - probability of a second yellow
# =============================================================================
cat("\n========== Result 1: second yellow probability ==========\n")

cat("\n[All bookings]\n")
booked |>
  summarise(n = n(), dismissals = sum(got_second),
            rate_pct = round(100*mean(got_second), 3)) |>
  as.data.frame() |> print()

# Restrict to bookings where a substitution decision is genuinely available.
analysis <- booked |> filter(mins_left >= MIN_AFTER)

cat("\n[At least", MIN_AFTER, "minutes remaining]\n")
analysis |>
  summarise(n = n(), dismissals = sum(got_second),
            rate_pct = round(100*mean(got_second), 3)) |>
  as.data.frame() |> print()

bt <- binom.test(sum(analysis$got_second), nrow(analysis))
cat("\n95% CI: [", round(100*bt$conf.int[1], 2), ",",
    round(100*bt$conf.int[2], 2), "] %\n")


# =============================================================================
# 4. Team-level panel in 5-minute bins
#
#    Unit of observation: one team, one match, one 5-minute bin. Each bin
#    carries the team's xG for minus its xG against, and a disciplinary state:
#
#      0 = no card        1 = carrying a booking      2 = a man down
#
#    The 5-minute resolution follows Badiella et al. (2023), who use the same
#    binning on the same five leagues, so that our estimates are directly
#    comparable to theirs.
# =============================================================================
cat("\n========== 4. Team-level 5-minute panel ==========\n")

teams_in_match <- events |>
  filter(teamId > 0) |>
  distinct(matchId, teamId) |>
  arrange(matchId, teamId)

# Map each team to its opponent, so xG against can be looked up.
opp_map <- teams_in_match |>
  group_by(matchId) |>
  filter(n() == 2) |>
  mutate(opp = rev(teamId)) |>
  ungroup()

team_book <- events |>
  filter(is_yellow, teamId > 0) |>
  group_by(matchId, teamId) |>
  summarise(t_book = min(match_time), .groups = "drop")

# State 2 covers both routes to playing a man down: a second yellow and a
# straight red. The cost of being short-handed does not depend on which.
team_dis <- events |>
  filter((is_second_yellow | is_red), teamId > 0) |>
  group_by(matchId, teamId) |>
  summarise(t_dis = min(match_time), .groups = "drop")

shot_bins <- shots |>
  mutate(bin = floor(match_time / 300)) |>
  group_by(matchId, teamId, bin) |>
  summarise(xg = sum(xG), .groups = "drop")

max_bin <- ceiling(max(match_end$match_end) / 300)

panel <- opp_map |>
  mutate(match_end = match_end$match_end[match(matchId, match_end$matchId)]) |>
  slice(rep(1:n(), each = max_bin)) |>
  group_by(matchId, teamId) |>
  mutate(bin = 0:(max_bin - 1)) |>
  ungroup() |>
  mutate(bin_mid = bin*300 + 150) |>
  filter(bin_mid < match_end) |>           # drop bins past the final whistle
  mutate(
    k_self  = paste(matchId, teamId, bin),
    k_opp   = paste(matchId, opp,    bin),
    k_team  = paste(matchId, teamId),
    xg_for  = shot_bins$xg[match(k_self, paste(shot_bins$matchId, shot_bins$teamId, shot_bins$bin))],
    xg_ag   = shot_bins$xg[match(k_opp,  paste(shot_bins$matchId, shot_bins$teamId, shot_bins$bin))],
    t_book  = team_book$t_book[match(k_team, paste(team_book$matchId, team_book$teamId))],
    t_dis   = team_dis$t_dis[ match(k_team, paste(team_dis$matchId,  team_dis$teamId))]
  ) |>
  mutate(
    # A bin with no shots is a genuine zero, not a missing value.
    xg_for  = ifelse(is.na(xg_for), 0, xg_for),
    xg_ag   = ifelse(is.na(xg_ag),  0, xg_ag),
    xg_diff = xg_for - xg_ag,
    minute  = bin_mid / 60,
    state   = case_when(
      !is.na(t_dis)  & bin_mid >= t_dis  ~ 2,
      !is.na(t_book) & bin_mid >= t_book ~ 1,
      TRUE                               ~ 0
    )
  )

cat("Bins:", nrow(panel), " / team-matches:", length(unique(panel$k_team)), "\n\n")
panel |>
  group_by(state) |>
  summarise(n = n(), mean_xg_diff = round(mean(xg_diff), 5), .groups = "drop") |>
  as.data.frame() |> print()
cat("\n(Validation gate: state 2 should have the lowest mean.)\n")


# =============================================================================
# 5. Result 2 - three-state regression
#
#    Model 1 controls only for match minute. Model 2 adds team-match fixed
#    effects by demeaning within team-match.
#
#    The fixed effects are not a refinement, they are the identification. Cards
#    are issued disproportionately to teams that are already defending and
#    under pressure, so a raw comparison of booked to unbooked periods
#    attributes the consequences of being outplayed to the card. Model 1 shows
#    what that confounding produces; Model 2 removes the team-match component
#    of it. The booked-state coefficient falls by roughly three quarters
#    between the two, and this is the single most important diagnostic in the
#    paper.
#
#    Coefficients are multiplied by 18 to convert from "per 5-minute bin" to
#    "per 90 minutes" (90 / 5 = 18).
# =============================================================================
cat("\n========== Result 2: three-state regression ==========\n")

panel <- panel |>
  mutate(s1 = as.numeric(state == 1),
         s2 = as.numeric(state == 2))

m1 <- lm(xg_diff ~ minute + s1 + s2, data = panel)
cat("\n[Model 1] Match minute only\n")
print(round(summary(m1)$coefficients, 6))

panel_w <- panel |>
  group_by(k_team) |>
  mutate(across(c(xg_diff, s1, s2, minute), ~ .x - mean(.x))) |>
  ungroup()

m2 <- lm(xg_diff ~ 0 + minute + s1 + s2, data = panel_w)
cat("\n[Model 2] Team-match fixed effects\n")
print(round(summary(m2)$coefficients, 6))

b    <- coef(m2)["s1"] * 18     # cost of carrying a booking, per 90
C    <- coef(m2)["s2"] * 18     # cost of being a man down, per 90
p_s1 <- summary(m2)$coefficients["s1", "Pr(>|t|)"]
p_s2 <- summary(m2)$coefficients["s2", "Pr(>|t|)"]

cat("\n--- Converted to per 90 minutes ---\n")
cat("Carrying a booking:", round(b, 4), " xG difference / 90   p =", signif(p_s1, 4), "\n")
cat("A man down        :", round(C, 4), " xG difference / 90   p =", signif(p_s2, 4), "\n")

# ---- Clustered standard errors ----------------------------------------------
#      Five-minute intervals within a team-match are correlated, and the two
#      teams in a match are mechanically related: xG difference is zero-sum
#      between them, so one team's series is the negative of the other's.
#      Team-match fixed effects absorb the team-match mean but not this
#      structure. We therefore cluster by match.
#
#      The coefficients are identical to m1 and m2 above; only the standard
#      errors change, and they change little, which is itself informative.
# =============================================================================

m1_cl <- feols(xg_diff ~ minute + s1 + s2,          data = panel, cluster = ~matchId)
m2_cl <- feols(xg_diff ~ minute + s1 + s2 | k_team, data = panel, cluster = ~matchId)

cat("\n[Model 1, clustered by match]\n"); print(summary(m1_cl))
cat("\n[Model 2, clustered by match]\n"); print(summary(m2_cl))

cat("\n--- Model 2, per 90 minutes, clustered ---\n")
for (v in c("s1", "s2")) {
  est <- coef(m2_cl)[v] * 18
  se  <- se(m2_cl)[v]   * 18
  cat(sprintf("%-3s %8.4f   SE %.4f   CI [%.3f, %.3f]   p = %.4g\n",
              v, est, se, est - 1.96*se, est + 1.96*se, pvalue(m2_cl)[v]))
}

b <- coef(m2_cl)["s1"] * 18
C <- coef(m2_cl)["s2"] * 18
p_s1 <- pvalue(m2_cl)["s1"]
p_s2 <- pvalue(m2_cl)["s2"]


# =============================================================================
# 6. Result 3 - expected cost and the substitution threshold
#
#    Leaving a booked player on the pitch carries two costs:
#      (i)  he plays the rest of the match in the booked state
#      (ii) with probability P he is dismissed, and the team then plays the
#           remainder a man down
#
#    Expected cost = b * (E_left / 90) + P * (C - b) * (E_down / 90)
#
#    The second term uses (C - b), not C. Both coefficients are measured
#    against the no-card state, and the first term already charges b over the
#    entire remaining window, including any minutes played after a dismissal.
#    Charging C over those same minutes as well would count them twice. What a
#    dismissal adds is the difference between the two states, not the level of
#    the man-down state.
#
#    Substituting avoids that cost but replaces the starter with someone of
#    lower value. Setting the two equal gives the threshold: the substitution
#    is justified only if the starter-to-replacement gap is smaller than
#
#      G*  =  | b + P * (C - b) * (E_down / E_left) |
#
#    Note that this is algebraically identical to dividing the total expected
#    cost by the exposure window. An earlier presentation of the same number
#    used the division form, which obscured the fact that the booked-state term
#    dominates. The expanded form above is preferred for that reason.
# =============================================================================
cat("\n========== Result 3: expected cost ==========\n")

p_dis  <- mean(analysis$got_second)
E_left <- mean(analysis$mins_left)
E_down <- mean((analysis$mins_left - analysis$expo_mins)[analysis$got_second == 1])

cost_booked <- b * (E_left / 90)
cost_dismis <- p_dis * (C - b) * (E_down / 90)
total_cost  <- cost_booked + cost_dismis

cat("\nP(dismissal | booked, >=", MIN_AFTER, "min remaining):", round(p_dis, 4), "\n")
cat("Mean minutes remaining after the card:", round(E_left, 1), "\n")
cat("Mean minutes played a man down, when dismissed:", round(E_down, 1), "\n")

cat("\nCost of the booked state:", round(cost_booked, 4), " xG\n")
cat("Cost via the dismissal path:", round(cost_dismis, 4), " xG\n")
cat("Total:", round(total_cost, 4), " xG\n")
cat("Dismissal share of total:", round(100*cost_dismis/total_cost, 1), "%\n")

cat("\n--- Threshold by minutes remaining:  G*(M) = b + P(M) x (C - b) x E_down(M)/M ---\n")
for (M in c(15, 30, 45, 60, 75)) {
  sub <- analysis |> filter(mins_left >= M - 7.5, mins_left < M + 7.5)
  if (nrow(sub) < 30) next
  p_M <- mean(sub$got_second)
  d_M <- mean((sub$mins_left - sub$expo_mins)[sub$got_second == 1])
  if (is.nan(d_M)) { d_M <- 0 }
  G <- b + p_M * (C - b) * (d_M / M)
  cat(sprintf("  %2d min left:  P = %.4f   G* = %+.4f xG/90   (n = %d)\n",
              M, p_M, G, nrow(sub)))
}
# Caution: P(M) is not monotonic across these strata. Each cell contains only
# a handful of dismissals, so the variation is sampling noise rather than a
# time profile, and the table should be read as a sensitivity check, not as an
# estimated function of M.


# =============================================================================
# 6.5 Robustness checks
# =============================================================================

# ---- (a) Confidence intervals on the per-90 coefficients ---------------------
cat("\n========== 6.5a Confidence intervals (per 90) ==========\n")
ci <- cbind("2.5 %"  = (coef(m2_cl) - 1.96 * se(m2_cl)) * 18,
            "97.5 %" = (coef(m2_cl) + 1.96 * se(m2_cl)) * 18)
print(round(ci, 4))

# ---- (b) Superseded ---------------------------------------------------------
#      An earlier version corrected the residual degrees of freedom by hand,
#      because demeaning does not tell lm() that 760 fixed effects had been
#      absorbed. Clustering by match (Section 5) subsumes that correction, and
#      also handles the dependence between the two teams in a match, which the
#      degrees-of-freedom adjustment did not address.

# ---- (c) Within-dismissal-match comparison ----------------------------------
#      Restricting to team-matches in which a dismissal actually occurred, with
#      team-match fixed effects, compares each team against itself before and
#      after going a man down. Without the fixed effects this block compared
#      different teams within the restricted sample, which is not the
#      comparison the text describes.
cat("\n========== 6.5c Matches containing a dismissal only ==========\n")
dis_teams <- panel |> filter(state == 2) |> distinct(k_team) |> pull(k_team)
sub_panel <- panel |> filter(k_team %in% dis_teams)
cat("Team-matches with a dismissal:", length(dis_teams), "\n\n")

sub_panel |> group_by(state) |>
  summarise(n = n(), mean_xg_diff = round(mean(xg_diff), 5), .groups = "drop") |>
  as.data.frame() |> print()

m3 <- feols(xg_diff ~ minute + s1 + s2 | k_team,
            data = sub_panel, cluster = ~matchId)
cat("\n[Model 3] Within dismissal matches, team-match FE, clustered by match\n")
print(summary(m3))
cat("\nA man down (per 90):", round(coef(m3)["s2"]*18, 4),
    "  p =", signif(pvalue(m3)["s2"], 4), "\n")

# =============================================================================
# 6.6 Starter-to-replacement value gap, first pass
#
#    The threshold above is only interpretable if we know how large a realistic
#    starter-to-replacement gap actually is. We estimate it with the action
#    value model: each player's per-90 count of each action type is weighted by
#    the corresponding regression coefficient and summed.
#
#    WARNING. The playing-time estimate used in this block is biased and the
#    result below should not be quoted. Minutes are approximated by the span
#    between a player's first and last recorded event, which understates time
#    on the pitch for substitutes, whose first touch may come well after they
#    enter. Because the understatement is proportionally largest for short
#    appearances, per-90 rates are inflated for exactly the players we classify
#    as replacements. Section 6.7 repeats the calculation without this bias and
#    is the version reported in the paper. This block is retained so that the
#    size and direction of the bias can be seen.
# =============================================================================
cat("\n========== 6.6 CDM starter vs replacement (biased minutes) ==========\n")

pm <- events |>
  filter(playerId > 0) |>
  group_by(playerId, matchId) |>
  summarise(mins = (max(match_time) - min(match_time))/60, .groups = "drop")

season_mins <- pm |>
  group_by(playerId) |>
  summarise(season_mins = sum(mins), n_matches = n(), .groups = "drop")

# Identify central defensive midfielders: players listed as midfielders whose
# mean action location falls in the deepest third of that group. This is a
# positional proxy, not a role label from the data provider.
avg_x <- events |>
  filter(playerId > 0, !is.na(start_x)) |>
  group_by(playerId) |>
  summarise(avg_x = mean(start_x), .groups = "drop")

mids <- avg_x |>
  mutate(role = players$role[match(playerId, players$playerId)]) |>
  filter(role == "Midfielder")

cut_x   <- quantile(mids$avg_x, 1/3, na.rm = TRUE)
cdm_ids <- mids |> filter(avg_x <= cut_x) |> pull(playerId)
cat("CDM candidates:", length(cdm_ids), "\n")

# Action value coefficients, as estimated in the original project. Each is the
# change in team xG difference per 90 associated with one additional action of
# that type.
BETA <- c(
  basic_pass      =  0.0034,  smart_pass      =  0.0968,
  key_pass        =  0.1180,  dribble_won     = -0.0048,
  cross_acc       =  0.0748,  shot            =  0.1760,
  tackle_won      = -0.0066,  tackle_fail     = -0.0189,
  aerial_won      =  0.0243,  aerial_fail     =  0.0210,
  turnover_own    = -0.0494,  turnover_opp    = -0.0362,
  foul_box        = -0.3500,  foul_def        =  0.0053
)

# Renamed to avoid name shadowing: assigning a column called season_mins inside
# mutate() would mask the lookup table of the same name.
mins_tbl <- season_mins |> rename(pid = playerId)

count_actions <- function(df) {
  df |> summarise(
    basic_pass   = sum(subEventName == "Simple pass", na.rm = TRUE),
    smart_pass   = sum(subEventName == "Smart pass",  na.rm = TRUE),
    key_pass     = sum(is_key_pass == 1,              na.rm = TRUE),
    dribble_won  = sum(subEventName == "Ground attacking duel" & duel_won, na.rm = TRUE),
    cross_acc    = sum(subEventName == "Cross" & is_accurate == 1,         na.rm = TRUE),
    shot         = sum(eventName == "Shot",           na.rm = TRUE),
    tackle_won   = sum(subEventName == "Ground defending duel" &  duel_won, na.rm = TRUE),
    tackle_fail  = sum(subEventName == "Ground defending duel" & !duel_won, na.rm = TRUE),
    aerial_won   = sum(subEventName == "Air duel" &  duel_won, na.rm = TRUE),
    aerial_fail  = sum(subEventName == "Air duel" & !duel_won, na.rm = TRUE),
    turnover_own = sum(is_dangerous_loss == 1 & start_x <  50, na.rm = TRUE),
    turnover_opp = sum(is_dangerous_loss == 1 & start_x >= 50, na.rm = TRUE),
    foul_box     = sum(eventName == "Foul" & start_x <  16, na.rm = TRUE),
    foul_def     = sum(eventName == "Foul" & start_x >= 16 & start_x < 50, na.rm = TRUE),
    n_fm = n_distinct(matchId), .groups = "drop")
}
# Duel outcomes use the duel_won column. Do not substitute tag 1801
# ("accurate") for this: it is a different field, and the natural sanity check
# on duels -- that win rates average 0.50 -- cannot detect the error, because
# duels are symmetric and the wrong field is symmetric too.

cdm_actions <- events |>
  filter(playerId %in% cdm_ids) |>
  group_by(playerId) |>
  count_actions() |>
  left_join(mins_tbl, by = c("playerId" = "pid")) |>
  filter(!is.na(season_mins), season_mins >= 270)   # at least three matches' worth

act_cols <- names(BETA)
val <- cdm_actions
for (a in act_cols) val[[paste0(a, "_p90")]] <- val[[a]] / val$season_mins * 90

val$value_p90 <- 0
for (a in act_cols) val$value_p90 <- val$value_p90 + BETA[a] * val[[paste0(a, "_p90")]]

val <- val |>
  mutate(player = players$shortName[match(playerId, players$playerId)])

cat("CDMs analysed:", nrow(val), "(270+ minutes)\n\n")
cat("Distribution of value per 90:\n")
print(round(summary(val$value_p90), 4))

cat("\n--- Classified by playing time ---\n")
starters <- val |> filter(season_mins >= 1800)   # 20 matches' worth or more
backups  <- val |> filter(season_mins <  900)    # fewer than 10 matches' worth

cat("Starters (1800+ min):", nrow(starters), " mean value per 90 =",
    round(mean(starters$value_p90), 4), "\n")
cat("Replacements (<900) :", nrow(backups),  " mean value per 90 =",
    round(mean(backups$value_p90), 4), "\n")

gap <- mean(starters$value_p90) - mean(backups$value_p90)
cat("\nGap (starter minus replacement):", round(gap, 4), "xG/90\n")

if (nrow(starters) > 1 & nrow(backups) > 1) {
  tt <- t.test(starters$value_p90, backups$value_p90)
  cat("t-test: t =", round(tt$statistic, 3),
      " p =", signif(tt$p.value, 4),
      " 95% CI [", round(tt$conf.int[1], 4), ",",
      round(tt$conf.int[2], 4), "]\n")
}

# Within-team comparison. This is more informative than the pooled comparison
# because value per 90 here is a plus-minus quantity: a replacement at a strong
# club can outscore a starter at a weak one without being the better player.
# Comparing within a squad holds team strength fixed.
cat("\n--- Within-team comparison ---\n")
val_team <- val |>
  # Team is taken from the player's first recorded event; players who moved
  # mid-season are assigned to one club only.
  mutate(teamId = events$teamId[match(playerId, events$playerId)]) |>
  group_by(teamId) |>
  filter(n() >= 2) |>
  mutate(rank_mins = rank(-season_mins, ties.method = "first"),
         is_top    = as.numeric(rank_mins == 1)) |>
  ungroup()

within_gap <- val_team |>
  group_by(teamId) |>
  summarise(top  = value_p90[is_top == 1][1],
            rest = mean(value_p90[is_top == 0]),
            .groups = "drop") |>
  mutate(gap = top - rest) |>
  filter(!is.na(gap))

cat("Teams compared:", nrow(within_gap), "\n")
cat("Mean within-team gap:", round(mean(within_gap$gap), 4), "xG/90\n")
cat("Median:", round(median(within_gap$gap), 4), "\n")
if (nrow(within_gap) > 1) {
  tt2 <- t.test(within_gap$gap)
  cat("t-test: p =", signif(tt2$p.value, 4),
      " 95% CI [", round(tt2$conf.int[1], 4), ",",
      round(tt2$conf.int[2], 4), "]\n")
}

threshold <- abs(total_cost / (E_left / 90))
cat("\n========================================\n")
cat("Threshold justifying substitution:", round(threshold, 4), "xG/90\n")
cat("Observed starter-replacement gap :", round(gap, 4), "xG/90  [BIASED - see 6.7]\n")
cat("Within-team gap                  :", round(mean(within_gap$gap), 4), "xG/90\n")
cat("========================================\n")


# =============================================================================
# 6.7 Starter-to-replacement gap, corrected
#
#    The bias in 6.6 is removed by using only appearances in which the player
#    demonstrably played the whole match, identified as a first event within
#    the opening four minutes and a last event after the 85th. For those
#    appearances the denominator is known to be 90 minutes, so no estimate of
#    playing time is required at all.
#
#    Classification into starter and replacement now uses the count of
#    appearances, which is unbiased, rather than estimated minutes.
#
#    One residual selection problem remains and should be stated in the paper:
#    requiring three full matches excludes players who never start, so the
#    surviving replacements are the better ones. The gap reported here is
#    therefore a lower bound.
# =============================================================================
cat("\n========== 6.7 Corrected: full-match appearances only ==========\n")

appear <- events |>
  filter(playerId > 0) |>
  group_by(playerId, matchId) |>
  summarise(first_t = min(match_time), last_t = max(match_time), .groups = "drop") |>
  mutate(full_match = first_t <= 240 & last_t >= 85*60)

n_app <- appear |> group_by(playerId) |>
  summarise(n_app = n(), n_full = sum(full_match), .groups = "drop")

fm_key <- appear |> filter(full_match) |> mutate(k = paste(playerId, matchId)) |> pull(k)

cdm_full <- events |>
  filter(playerId %in% cdm_ids, paste(playerId, matchId) %in% fm_key) |>
  group_by(playerId) |>
  count_actions() |>
  left_join(n_app, by = "playerId") |>
  filter(n_fm >= 3)

# Each retained appearance is exactly 90 minutes, so the per-90 rate is simply
# the total count divided by the number of full matches.
cdm_full$value_p90 <- 0
for (a in names(BETA))
  cdm_full$value_p90 <- cdm_full$value_p90 + BETA[a] * (cdm_full[[a]] / cdm_full$n_fm)

cdm_full <- cdm_full |>
  mutate(player = players$shortName[match(playerId, players$playerId)],
         grp = ifelse(n_app >= 25, "starter",
                      ifelse(n_app <= 15, "replacement", "intermediate")))

cat("CDMs analysed:", nrow(cdm_full), "(3+ full matches)\n\n")
cdm_full |> group_by(grp) |>
  summarise(n = n(),
            mean_appearances = round(mean(n_app), 1),
            mean_full_matches = round(mean(n_fm), 1),
            value_p90 = round(mean(value_p90), 4), .groups = "drop") |>
  as.data.frame() |> print()

s_grp <- cdm_full |> filter(grp == "starter")
r_grp <- cdm_full |> filter(grp == "replacement")
if (nrow(s_grp) > 1 & nrow(r_grp) > 1) {
  g2 <- mean(s_grp$value_p90) - mean(r_grp$value_p90)
  t2 <- t.test(s_grp$value_p90, r_grp$value_p90)
  cat("\nCorrected gap:", round(g2, 4), " xG/90   p =", signif(t2$p.value, 4),
      "  95% CI [", round(t2$conf.int[1], 4), ",", round(t2$conf.int[2], 4), "]\n")
  cat("Compared with threshold", round(threshold, 4), "->",
      ifelse(g2 > threshold, "gap exceeds threshold", "gap below threshold"), "\n")
}
# Note for the write-up: the three groups are not ordered by playing time.
# Value per 90 does not rank starters above replacements in this sample, which
# suggests the action value model is not measuring the attributes on which
# managers select -- positioning, off-ball defending, tactical discipline and
# the severity rather than the count of errors are all outside it.


# =============================================================================
# 6.8 Sensitivity of the threshold
#
#    The threshold depends on b, whose p-value is 0.487 once standard errors
#    are clustered by match. Since b is not distinguishable from zero, the
#    threshold is not pinned down either, and reporting a single value would
#    overstate what the data support.
# =============================================================================
cat("\n========== 6.8 Threshold sensitivity ==========\n")
thr <- function(bb, CC) abs(bb*(E_left/90) + p_dis*(CC - bb)*(E_down/90)) / (E_left/90)
cat("Point estimate (b =", round(b, 3), "):", round(thr(b, C), 4), "\n")
cat("Assuming b = 0 (no booked-state effect):", round(thr(0, C), 4), "\n")
cat("Lower bounds of b and C:", round(thr(ci["s1","2.5 %"], ci["s2","2.5 %"]), 4), "\n")

# How much data would settle it? To reach p < 0.05 at the observed effect size,
# the clustered standard error on b must fall to |b| / 1.96.
se_b       <- se(m2_cl)["s1"] * 18
se_needed  <- abs(b) / 1.96
n_factor   <- (se_b / se_needed)^2
cat("\nStandard error on b:", round(se_b, 4),
    " / required for p < 0.05:", round(se_needed, 4), "\n")
cat("Sample size multiple required:", round(n_factor, 1),
    " -> approximately", round(380 * n_factor), "matches\n")
cat("(Five leagues over two seasons would reach this.)\n")

# =============================================================================
# 6.9 Second yellow probability for central defensive midfielders only
#
#    The expected cost in Section 6 applies a league-wide dismissal rate to a
#    position-specific decision. CDMs contest more duels than the average
#    player, so their risk need not match the league rate. This block
#    recomputes P on the CDM subset alone and rebuilds the threshold from it.
#
#    The subset is small, so the interval is wide. This is reported as a
#    descriptive comparison, not as a test of whether CDMs differ.
#
#    cdm_ids is defined in Section 6.6 and reused here unchanged.
# =============================================================================
cat("\n========== 6.9 Second yellow probability, CDMs only ==========\n")

analysis_cdm <- analysis |> filter(playerId %in% cdm_ids)
analysis_oth <- analysis |> filter(!(playerId %in% cdm_ids))

cat("CDM booked player-matches (>=", MIN_AFTER, "min remaining):",
    nrow(analysis_cdm), "\n")
cat("Dismissals:", sum(analysis_cdm$got_second), "\n")

p_dis_cdm <- mean(analysis_cdm$got_second)
bt_cdm    <- binom.test(sum(analysis_cdm$got_second), nrow(analysis_cdm))

cat("Rate:", round(100*p_dis_cdm, 3), "%   95% CI: [",
    round(100*bt_cdm$conf.int[1], 2), ",",
    round(100*bt_cdm$conf.int[2], 2), "] %\n")
cat("League-wide rate for comparison:", round(100*p_dis, 3), "%\n")

# ---- Is the CDM rate distinguishable from everyone else's? -------------------
tab <- matrix(
  c(sum(analysis_cdm$got_second), nrow(analysis_cdm) - sum(analysis_cdm$got_second),
    sum(analysis_oth$got_second), nrow(analysis_oth) - sum(analysis_oth$got_second)),
  nrow = 2, byrow = TRUE,
  dimnames = list(c("CDM", "Other"), c("dismissed", "not")))
cat("\n")
print(tab)
cat("Fisher exact p =", signif(fisher.test(tab)$p.value, 4), "\n")

# ---- Expected cost and threshold, rebuilt with the CDM-specific rate ---------
E_left_cdm <- mean(analysis_cdm$mins_left)
E_down_cdm <- mean((analysis_cdm$mins_left -
                      analysis_cdm$expo_mins)[analysis_cdm$got_second == 1])
if (is.nan(E_down_cdm)) E_down_cdm <- E_down   # fallback if no CDM dismissals

cost_booked_cdm <- b * (E_left_cdm / 90)
cost_dismis_cdm <- p_dis_cdm * (C - b) * (E_down_cdm / 90)
total_cdm       <- cost_booked_cdm + cost_dismis_cdm
thr_cdm         <- abs(total_cdm / (E_left_cdm / 90))

cat("\n--- Recomputed with the CDM-specific rate ---\n")
cat("E_left:", round(E_left_cdm, 1), "  E_down:", round(E_down_cdm, 1), "\n")
cat("Booked state:", round(cost_booked_cdm, 4),
    "  Dismissal path:", round(cost_dismis_cdm, 4), "\n")
cat("Total:", round(total_cdm, 4),
    "  dismissal share:", round(100*cost_dismis_cdm/total_cdm, 1), "%\n")
cat("Threshold G*:", round(thr_cdm, 4),
    "  (league-wide version:", round(threshold, 4), ")\n")


# =============================================================================
# 7. Figure 8
#    Grey and gold share one time slope with a booked dummy, mirroring the
#    player-level specification: both states are eleven against eleven, so a
#    common trend is a defensible restriction. The man-down series is fitted on
#    its own data, because a team playing ten against eleven is not the same
#    game and the eleven-a-side trend should not be imposed on it.
# =============================================================================
m01 <- lm(xg_diff ~ minute + s1, data = filter(panel, state != 2, minute <= 92))
m2  <- lm(xg_diff ~ minute,      data = filter(panel, state == 2, minute <= 92))

cat("\n--- Figure 8 fitted lines ---\n")
cat("States 0/1 common slope :", round(coef(m01)["minute"], 6), "\n")
cat("Booked offset           :", round(coef(m01)["s1"], 6),
    " (", round(coef(m01)["s1"]*18, 4), " per 90 )\n")
cat("Man-down own slope      :", round(coef(m2)["minute"], 6), "\n")

rng <- function(st) range(panel$minute[panel$state == st])
g01 <- bind_rows(
  data.frame(state = 0, minute = seq(rng(0)[1], rng(0)[2], length.out = 200), s1 = 0),
  data.frame(state = 1, minute = seq(rng(1)[1], rng(1)[2], length.out = 200), s1 = 1))
g01$fit <- predict(m01, newdata = g01)

g2 <- data.frame(state = 2, minute = seq(rng(2)[1], rng(2)[2], length.out = 200))
g2$fit <- predict(m2, newdata = g2)

fit_lines <- bind_rows(g01[, c("state","minute","fit")], g2[, c("state","minute","fit")])
fit_lines$state_lab <- factor(fit_lines$state, levels = c(0,1,2),
                              labels = c("No card","After first yellow","A man down"))

pts <- panel |>
  filter(minute <= 92) |>
  mutate(bin5 = floor(minute/5)*5,
         state_lab = factor(state, levels = c(0,1,2),
                            labels = c("No card","After first yellow","A man down"))) |>
  group_by(state_lab, bin5) |>
  summarise(m = mean(xg_diff), n = n(), .groups = "drop") |>
  filter(n >= MIN_BIN_N)

gx    <- 82
g_hi  <- coef(m01)[1] + coef(m01)["minute"]*gx
g_lo  <- g_hi + coef(m01)["s1"]
g_lab <- sprintf("uncontrolled gap\n%+.3f xG/90", coef(m01)["s1"]*18)

p <- ggplot() +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_point(data = pts, aes(bin5, m, colour = state_lab), size = 2.6) +
  geom_line(data = fit_lines, aes(minute, fit, colour = state_lab), linewidth = 1.1) +
  annotate("segment", x = gx, xend = gx, y = g_lo, yend = g_hi,
           colour = "firebrick4", linewidth = 0.6,
           arrow = arrow(ends = "both", length = unit(0.16, "cm"))) +
  annotate("text", x = gx + 1.5, y = (g_lo + g_hi)/2, label = g_lab,
           hjust = 0, size = 3.1, colour = "firebrick4", lineheight = 0.95) +
  scale_colour_manual(values = c("No card"            = "grey45",
                                 "After first yellow" = "goldenrod3",
                                 "A man down"         = "firebrick")) +
  scale_x_continuous(limits = c(-2, 106), breaks = seq(0, 90, 15)) +
  labs(x = "match minute (5-minute bins)",
       y = "Team xG difference per 5-minute bin",
       colour = NULL,
       title = "Team xG difference by disciplinary state",
       subtitle = "La Liga 2017-18") +
  theme_minimal(base_size = 13) +
  theme(legend.position = "top")

print(p)
ggsave(file.path(DATA_DIR, "fig_three_state.png"), p,
       width = 9, height = 5, dpi = 300, bg = "white")

# =============================================================================
# Session information, for reproducibility
# =============================================================================
cat("\n========== Session information ==========\n")
print(sessionInfo())