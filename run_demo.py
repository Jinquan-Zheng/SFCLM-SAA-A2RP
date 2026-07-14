# Demo script: run SFCLM+A2RP on bundled sample scenarios
from pathlib import Path
import os
os.chdir(Path(__file__).resolve().parent)

# SFCLM core functions: data loading, two-stage stochastic model, and A2RP
import ast
import os
import random
from collections import defaultdict
from pathlib import Path

import gurobipy as gp
import numpy as np
import pandas as pd
from gurobipy import GRB
from scipy.sparse import lil_matrix
from scipy.stats import t
from tqdm import tqdm

# Reproducibility seed used by A2RP scenario splits
RANDOM_SEED = 20251119


################### Data loading and preprocessing #################

def load_scenarios_from_folder(folder_path):
    """Load scenario CSV paths from a folder (exclude SAA manifest.csv)."""
    scenario_files = [
        os.path.join(folder_path, file)
        for file in os.listdir(folder_path)
        if file.endswith(".csv") and file != "manifest.csv"
    ]
    print(f"Number of scenario files loaded: {len(scenario_files)}")
    return scenario_files


def load_shortest_paths_from_file(od_file):
    """Load shortest-path node sequences for trips with final_mode == Multi-hop."""
    od_df = pd.read_csv(od_file)
    # Keep only trips that chose the multi-hop ride-pooling alternative
    od_df = od_df[od_df["final_mode"] == "Multi-hop"]
    short_paths = []

    for path in od_df["shortest_path"]:
        if pd.isna(path):
            continue  # Skip missing path records
        # Parse path string safely into a list of node IDs
        short_paths.append(ast.literal_eval(path))

    return short_paths


def generate_edge_flows(short_paths):
    """Convert paths into undirected edge flows as (start, end, flow)."""
    edge_flow = defaultdict(int)
    for path in short_paths:
        for i in range(len(path) - 1):
            # Ordered tuple so each undirected edge is unique
            edge = tuple(sorted((path[i], path[i + 1])))
            edge_flow[edge] += 1

    return [(e[0], e[1], flow) for e, flow in edge_flow.items()]


def get_filtered_candidate_stations(edges_and_flows, candidate_stations):
    """Return candidate stations that appear as endpoints of the given edges."""
    station_set = set()
    for start, end, _ in edges_and_flows:
        if start in candidate_stations:
            station_set.add(start)
        if end in candidate_stations:
            station_set.add(end)
    return station_set


def calculate_flow_for_solution(edges_and_flows, solution):
    """Compute total flow captured by a given location set."""
    total_flow = 0
    solution_set = set(solution)
    for start, end, flow in edges_and_flows:
        if start in solution_set or end in solution_set:
            total_flow += flow
    return total_flow


def create_F_and_A(edges_and_flows, candidate_stations):
    """Build objective vector F and sparse incidence matrix A for a single-scenario FCLM."""
    filtered_edges_and_flows = [
        (start, end, flow)
        for (start, end, flow) in edges_and_flows
        if start in candidate_stations or end in candidate_stations
    ]

    involved_stations = set()
    for start, end, _ in filtered_edges_and_flows:
        if start in candidate_stations:
            involved_stations.add(start)
        if end in candidate_stations:
            involved_stations.add(end)

    sorted_involved_stations = sorted(involved_stations)
    station_to_row = {station: idx for idx, station in enumerate(sorted_involved_stations)}

    num_stations = len(sorted_involved_stations)
    num_edges = len(filtered_edges_and_flows)

    A = lil_matrix((num_stations, num_edges), dtype=np.int8)
    F = np.zeros(num_edges)

    for j, (start, end, flow) in enumerate(filtered_edges_and_flows):
        F[j] = flow
        if start in station_to_row:
            A[station_to_row[start], j] = 1
        if end in station_to_row:
            A[station_to_row[end], j] = 1

    return F, A.tocsr(), sorted_involved_stations


def gurobi_progress_bar(model, where, progress):
    """Callback: advance progress bar whenever an integer solution is found."""
    if where == gp.GRB.Callback.MIPSOL:
        progress.update(1)


def flow_capture_model(F, A, p):
    """Solve a deterministic flow-capturing MIP (legacy helper; kept for completeness)."""
    model = gp.Model("FlowCapturing")
    model.setParam("OutputFlag", 0)
    x = model.addVars(A.shape[0], vtype=GRB.BINARY, name="x")
    y = model.addVars(len(F), vtype=GRB.BINARY, name="y")

    model.setObjective(gp.quicksum(F[j] * y[j] for j in range(len(F))), GRB.MAXIMIZE)

    A_csr = A.tocsr()
    for j in range(len(F)):
        rows = A_csr[:, j].nonzero()[0]
        model.addConstr(y[j] <= gp.quicksum(x[i] for i in rows), name=f"edge_{j}_capture")

    model.addConstr(gp.quicksum(x[i] for i in range(A.shape[0])) == p, name="facility_limit")

    with tqdm(total=1, desc="Gurobi Optimization Progress") as progress:
        model.optimize(lambda model, where: gurobi_progress_bar(model, where, progress))

    if model.status == GRB.OPTIMAL:
        return model.objVal, [i for i in range(A.shape[0]) if x[i].X > 0.5]
    return None, []


def preprocess_scenarios(scenario_files, candidate_stations):
    """Preprocess all scenarios into edge-flow lists and the active candidate subset."""
    scenario_data = []
    filtered_candidate_stations = set()
    for scenario_file in tqdm(scenario_files, desc="Preprocessing Scenarios"):
        short_paths = load_shortest_paths_from_file(scenario_file)
        edges_and_flows = generate_edge_flows(short_paths)
        scenario_data.append(edges_and_flows)

        stations_in_edges = get_filtered_candidate_stations(edges_and_flows, candidate_stations)
        filtered_candidate_stations.update(stations_in_edges)
    print(f"Number of candidate stations appearing in scenarios: {len(filtered_candidate_stations)}")
    return scenario_data, filtered_candidate_stations


############ Two-stage stochastic SFCLM (Gurobi) ###########################

def solve_two_stage_stochastic_model(filtered_candidate_stations, scenario_data, p):
    """Solve the two-stage SFCLM: choose p locations to maximize expected captured flow."""
    model = gp.Model("TwoStageStochasticFlowCapturing")
    model.setParam("OutputFlag", 0)

    candidate_stations_list = sorted(filtered_candidate_stations)
    num_stations = len(candidate_stations_list)
    station_to_index = {station: idx for idx, station in enumerate(candidate_stations_list)}

    # First-stage location variables x_i
    x_vars = model.addVars(num_stations, vtype=GRB.BINARY, name="x")

    num_scenarios = len(scenario_data)
    y_vars = {}
    objective_terms = []

    # Second-stage capture variables and expected-value objective
    for s in range(num_scenarios):
        edges_and_flows = scenario_data[s]
        y_vars_s = {}
        for start, end, flow in edges_and_flows:
            nodes = []
            if start in station_to_index:
                nodes.append(station_to_index[start])
            if end in station_to_index:
                nodes.append(station_to_index[end])
            if nodes:
                y_var = model.addVar(vtype=GRB.BINARY, name=f"y_{s}_{start}_{end}")
                y_vars_s[(start, end)] = (y_var, flow, nodes)
                # Average over scenarios to maximize expected captured flow
                objective_terms.append((1 / num_scenarios) * flow * y_var)
        y_vars[s] = y_vars_s

    model.setObjective(gp.quicksum(objective_terms), GRB.MAXIMIZE)
    model.addConstr(gp.quicksum(x_vars[i] for i in range(num_stations)) == p, name="FacilityLimit")

    # Flow can be captured only if at least one incident selected station exists
    for s in range(num_scenarios):
        for (start, end), (y_var, flow, nodes) in y_vars[s].items():
            model.addConstr(
                y_var <= gp.quicksum(x_vars[i] for i in nodes),
                name=f"FlowCapture_{s}_{start}_{end}",
            )

    model.optimize()

    STATUS_MAP = {
        1: "LOADED",
        2: "OPTIMAL",
        3: "INFEASIBLE",
        4: "INF_OR_UNBD",
        5: "UNBOUNDED",
        6: "CUTOFF",
        7: "ITERATION_LIMIT",
        8: "NODE_LIMIT",
        9: "TIME_LIMIT",
        10: "SOLUTION_LIMIT",
        11: "INTERRUPTED",
        12: "NUMERIC",
        13: "SUBOPTIMAL",
    }
    solve_stats = {
        "num_vars": model.NumVars,
        "num_constrs": model.NumConstrs,
        "solve_time": model.Runtime,
        "MIPGap": model.MIPGap if model.SolCount > 0 else None,
        "status": STATUS_MAP.get(model.status, str(model.status)),
    }

    if model.status == GRB.OPTIMAL:
        selected_facilities = [
            candidate_stations_list[i] for i in range(num_stations) if x_vars[i].X > 0.5
        ]
        total_captured_flow = model.ObjVal
    else:
        selected_facilities = None
        total_captured_flow = None
    return selected_facilities, total_captured_flow, solve_stats


def evaluate_solution(edges_and_flows, candidate_solution):
    """Evaluate captured flow of a candidate location set on one scenario."""
    return calculate_flow_for_solution(edges_and_flows, candidate_solution)


def solve_subproblem(num, indices, candidate_stations, scenario_data, p):
    """Solve the two-stage SFCLM on a scenario subset (used inside A2RP)."""
    sub_scenario_data = [scenario_data[i] for i in indices]
    selected_facilities, total_captured_flow, _ = solve_two_stage_stochastic_model(
        candidate_stations, sub_scenario_data, p
    )
    print(f"subproblem{num} Optimal solution stations: {selected_facilities}")
    print(f"subproblem{num} Captured flow: {total_captured_flow}")
    return total_captured_flow, selected_facilities


def a2rp_evaluation(best_solution, candidate_stations, scenario_data, L, p, seed=RANDOM_SEED):
    """Average two-replication procedure (A2RP) for optimality-gap assessment."""
    random.seed(seed)
    indices = list(range(L))
    random.shuffle(indices)
    I1, I2 = indices[: L // 2], indices[L // 2 :]

    z_star_I1, x_star_I1 = solve_subproblem(1, I1, candidate_stations, scenario_data, p)
    z_star_I2, x_star_I2 = solve_subproblem(2, I2, candidate_stations, scenario_data, p)

    z_bar_candidate_I1 = (1 / len(I1)) * sum(
        evaluate_solution(scenario_data[i], best_solution) for i in I1
    )
    print("optimal solution Captured subproblem1 flow:", round(z_bar_candidate_I1, 6))
    z_bar_candidate_I2 = (1 / len(I2)) * sum(
        evaluate_solution(scenario_data[i], best_solution) for i in I2
    )
    print("optimal solution Captured subproblem2 flow:", round(z_bar_candidate_I2, 6))

    G_I1 = z_star_I1 - z_bar_candidate_I1
    G_I2 = z_star_I2 - z_bar_candidate_I2
    print("G_I1:", round(G_I1, 6))
    print("G_I2:", round(G_I2, 6))

    s2_I1 = sum(
        (
            evaluate_solution(scenario_data[i], x_star_I1)
            - evaluate_solution(scenario_data[i], best_solution)
            - G_I1
        )
        ** 2
        for i in I1
    ) / (len(I1) - 1)

    s2_I2 = sum(
        (
            evaluate_solution(scenario_data[i], x_star_I2)
            - evaluate_solution(scenario_data[i], best_solution)
            - G_I2
        )
        ** 2
        for i in I2
    ) / (len(I2) - 1)

    print("s2_I1:", round(s2_I1, 6))
    print("s2_I2:", round(s2_I2, 6))

    G_I = (G_I1 + G_I2) / 2
    s2_I = (s2_I1 + s2_I2) / 2

    t_value = t.ppf(0.975, df=(L // 2) - 1)
    ci_upper = G_I + t_value * np.sqrt(s2_I / L)
    return round(G_I, 6), round(ci_upper, 6)


# Demo / main: solve SFCLM on sample scenarios and assess solution quality with A2RP
import csv
from pathlib import Path

# ===========================
# Paths relative to the repository root (parent of notebooks/)
REPO_ROOT = Path(__file__).resolve().parent
if REPO_ROOT.name == "notebooks":
    REPO_ROOT = REPO_ROOT.parent

# Demo defaults (small run for GitHub sample data).
# Paper-scale settings: p = 1..191 on full SAA folders total(100)..total(2000).
DEMO_MODE = True
P_MAX = 20 if DEMO_MODE else 191
SCENARIO_DIR = REPO_ROOT / "data" / "sample_scenarios"
CANDIDATE_STATIONS_FILE = REPO_ROOT / "data" / "candidate_stations_2019.csv"
OUTPUT_DIR = REPO_ROOT / "results"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# Relative optimality-gap tolerance used in the adaptive sample-size assessment
XI = 0.001

if __name__ == "__main__":
    candidate_stations_df = pd.read_csv(CANDIDATE_STATIONS_FILE, encoding="utf-8-sig")
    candidate_stations = set(candidate_stations_df["label"])

    print(f"Scenario directory: {SCENARIO_DIR}")
    print(f"Candidate stations file: {CANDIDATE_STATIONS_FILE}")
    print(f"Demo mode={DEMO_MODE}, p = 1..{P_MAX}, xi={XI}")

    scenario_files = load_scenarios_from_folder(str(SCENARIO_DIR))
    scenario_data, filtered_candidate_stations = preprocess_scenarios(
        scenario_files, candidate_stations
    )

    if not filtered_candidate_stations:
        raise RuntimeError(
            "No candidate stations found in the scenario data. Please check the input files."
        )

    max_effective_p = len(filtered_candidate_stations)
    remaining_candidates = sorted(candidate_stations - filtered_candidate_stations)
    print(
        f"Max effective p: {max_effective_p}; remaining unused candidates: {len(remaining_candidates)}"
    )

    sample_size = len(scenario_data)
    CI = []
    max_flow_solution = None
    max_flow_value = None
    max_flow_stats = None

    for p in range(1, P_MAX + 1):
        if p <= max_effective_p:
            selected_facilities, total_captured_flow, solve_stats = solve_two_stage_stochastic_model(
                filtered_candidate_stations, scenario_data, p
            )

            if selected_facilities is None:
                print(f"p={p}: solve failed; no optimal location set found.")
                continue

            if p == max_effective_p:
                max_flow_solution = selected_facilities
                max_flow_value = total_captured_flow
                max_flow_stats = solve_stats

            L = len(scenario_data)
            G_I, ci_upper = a2rp_evaluation(
                selected_facilities, candidate_stations, scenario_data, L, p
            )
        else:
            # Beyond the active candidate set, captured flow cannot increase further
            extra_needed = p - max_effective_p
            selected_facilities = max_flow_solution + remaining_candidates[:extra_needed]
            total_captured_flow = max_flow_value
            ci_upper = 0.0
            solve_stats = max_flow_stats

        ROG = ci_upper / total_captured_flow if total_captured_flow > 0 else 0.0
        print(f"p={p}: captured_flow={total_captured_flow:.6f}, CI_upper={ci_upper:.6f}, ROG={ROG:.6f}")

        CI.append(
            (
                p,
                str(selected_facilities),
                float(total_captured_flow),
                float(ci_upper),
                float(ROG),
                solve_stats["num_vars"],
                solve_stats["num_constrs"],
                solve_stats["MIPGap"],
                round(solve_stats["solve_time"], 4),
                solve_stats["status"],
            )
        )

    output_file = OUTPUT_DIR / f"{sample_size}_SFCLM_results(1-{P_MAX}).csv"
    with open(output_file, "w", newline="", encoding="utf-8") as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(
            [
                "p",
                "selected_facilities",
                "captured_flow",
                "CI_upper",
                "ROG",
                "num_vars",
                "num_constrs",
                "MIPGap",
                "solve_time_sec",
                "gurobi_status",
            ]
        )
        writer.writerows(CI)

    # Separator between per-p A2RP outputs and the final sample-size summary
    print("=" * 11)

    # Unified adaptive sample-size check: M is adequate only if ROG < xi for ALL p;
    # if any p has ROG >= xi, the current sample size is insufficient.
    inadequate_ps = [row[0] for row in CI if row[4] >= XI]
    if inadequate_ps:
        print(
            f"Sample-size check: M={sample_size} is NOT adequate "
            f"(ROG >= xi={XI} for p in {inadequate_ps}). Consider a larger SAA sample size."
        )
    else:
        print(
            f"Sample-size check: M={sample_size} is adequate "
            f"(ROG < xi={XI} for all tested p = 1..{P_MAX})."
        )

    print(f"Results for sample size M={sample_size} saved to {output_file}")

