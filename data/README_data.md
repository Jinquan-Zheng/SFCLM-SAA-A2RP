# Sample data notes

## Contents

- `candidate_stations_2019.csv`: 191 metro stations built in Shenzhen by the end of 2019 (candidate transfer locations). The `Name` column uses official English station names; key column `label` is the road-network node ID.
- `sample_scenarios/`: 15 scenario CSV files sampled from a representative weekday SAA set. Each file contains about 1,000 trip records.

## Scenario CSV schema

| Column          | Description                                                      |
| --------------- | ---------------------------------------------------------------- |
| origin          | Trip origin node ID                                              |
| destination     | Trip destination node ID                                         |
| shortest_path   | Solo-ride shortest-path node sequence (stringified list)         |
| travel_distance | Path distance                                                    |
| scenario_type   | Intra-day demand regime label (e.g., morning/afternoon/off-peak) |
| travel_time     | Path travel time                                                 |
| final_mode      | Simulated mode choice; SFCLM uses rows with `Multi-hop` only     |

These files are **demo inputs**. They are sufficient to run the notebook end-to-end, but numerical results will **not** match the full paper experiments (which use much larger SAA sample sizes).

## Full experimental data

The complete SAA scenario sets used in the paper are available from the corresponding author upon reasonable request.
