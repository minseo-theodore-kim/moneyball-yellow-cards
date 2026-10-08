# The Price of Playing Booked: analysis code

R code for the extension analyses in:

**The Price of Playing Booked: Yellow Cards and Central Defensive
Midfielder Value in La Liga**
Minseo (Theodore) Kim, Alexandre Meiler, Milo Lin, Yuvaan Pandey,
Harish Anand. Wharton Moneyball Academy, 2026.
Submitted to the Wharton Sports Analytics Journal.

---

## Scope

This repository covers **Sections 3.7 to 3.9** of the paper:

- the probability that a booked player receives a second yellow card
- the team-level three-state panel and the fixed-effects regression
- the expected cost of leaving a booked player on the pitch, and the
  substitution threshold derived from it
- the starter-to-replacement value gap
- the cross-league replication of the second-yellow rate

**It does not cover Sections 3.1 to 3.6.** The duel-behaviour tests and
the action value regression reported there were produced during the
competition stage of the project, in code that is not included here.

One consequence is worth stating. The `BETA` vector in Section 6.6 of
`Analysis_secondYellow.R` holds the fourteen action weights estimated by
that earlier regression. Those values are hard-coded here rather than
re-estimated, so this repository reproduces everything downstream of
them but not the weights themselves.

---

## Data

The analyses use the public Wyscout event data released by
Pappalardo et al. (2019):

> Pappalardo, L., Cintia, P., Rossi, A., Massucco, E., Ferragina, P.,
> Pedreschi, D., & Giannotti, F. (2019). A public data set of
> spatio-temporal match events in soccer competitions.
> *Scientific Data*, 6, 236. https://doi.org/10.1038/s41597-019-0247-7

The data files are not included in this repository. Download them and
place them in a `data/` folder at the repository root:

```
data/
├── events_Spain.csv     (required by both scripts)
├── players.csv          (required by Analysis_secondYellow.R)
├── events_Italy.csv     (required by part_3.9_research.R)
├── events_France.csv
├── events_Germany.csv
└── events_England.csv
```

`events_Spain.csv` is a preprocessed file carrying `is_yellow`,
`is_second_yellow`, `is_red` and `match_time` as columns. The other four
are the raw release files, which carry a `tags` string and
`eventSec`/`matchPeriod` instead. Both schemas are handled.

If your data lives elsewhere, change the `DATA_DIR` line at the top of
each script.

---

## Requirements

R 4.0 or later, with:

```r
install.packages(c("tidyverse", "fixest"))
```

---

## How to run

Run in this order:

```r
source("Analysis_secondYellow.R")   # Sections 3.7, 3.8
source("part_3.9_research.R")       # Section 3.9
```

Each script runs top to bottom and prints its results to the console.
Two figures are written to `DATA_DIR`:

| File | Paper |
|---|---|
| `fig_three_state.png` | Figure 8 |
| `fig_cross_league_sy.png` | Figure 9 |

---

## Validation gates

Both scripts stop or warn if the data do not reproduce known counts.
These are not cosmetic: an earlier version of the tag parser had the
yellow and red card identifiers reversed, and the gate in
`part_3.9_research.R` is what caught it.

`Analysis_secondYellow.R` checks, against the La Liga 2017-18 file:

```
first yellows      1,863
second yellows        42
straight reds         29
matches              380
predicted / actual goals ~ 1.00
state 2 has the lowest mean xG difference
```

`part_3.9_research.R` checks that the tag parser reproduces the stored
columns in `events_Spain.csv` exactly, and then that the full pipeline
reproduces Section 3.7:

```
1,365 bookings with at least 15 minutes remaining
   38 dismissals
2.784%
```

If the Spain row does not match, the script stops rather than printing
figures for the other four leagues.

---

## Headline results reproduced

```
Carrying a booking, team level, team-match FE   -0.064 per 90  (p = 0.487)
A man down, team level, team-match FE           -0.718 per 90  (p = 0.006)
P(second yellow | booked, 15+ min remaining)     2.78%  [1.98, 3.80]
Substitution threshold G*                        0.073 per 90
Starter-to-replacement gap                       0.015 per 90  (p = 0.816)
```

Standard errors are clustered by match throughout.

---

## Note on tooling

This repository contains only the two extension scripts, which the first author wrote after the Academy with coding assistance from Claude (Anthropic). The competition-stage code for the expected-goals and action value models is not included here (see the scope note above). All analyses were run by the authors in R, and every number reported in the paper comes from these scripts. The paper includes a fuller AI use statement.
