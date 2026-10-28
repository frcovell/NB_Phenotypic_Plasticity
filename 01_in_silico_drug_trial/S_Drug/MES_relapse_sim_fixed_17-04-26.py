import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
import itertools
import pandas as pd
import pickle
import os
from functools import partial
from scipy import integrate
import numpy as np
import sqlite3
import itertools
import matplotlib.pyplot as plt
import time as pytime


# =============================================================================
# PARAMETER SPACE GENERATION
# =============================================================================
# Creates every combination of the supplied parameters and saves:
#   - one CSV summary (InitialParamSpace.csv)
#   - one .pkl file per combination (Combo1.pkl, Combo2.pkl, ...)

def ParamSpace(pklSave, save_dir, n, m, s,
               K1_N, K1_I, K1_S,
               K2_N, K2_I, K2_S,
               K3_NS, K3_NI, K3_IN, K3_IS, K3_SI, K3_SN,
               Deg, Ind):

    keys = ['N', 'I', 'S',
            'K1N', 'K1I', 'K1S',
            'K2N', 'K2I', 'K2S',
            'K3NS', 'K3NI', 'K3IN', 'K3IS', 'K3SI', 'K3SN',
            'ki', 'kw']

    data = []

    # --- Iterate over all parameter combinations ---
    for j in range(len(n)):
        N, I, S = n[j], m[j], s[j]                         # Initial cell populations

        for k in range(len(K1_N)):
            K1N, K1I, K1S = K1_N[k], K1_I[k], K1_S[k]    # Intrinsic growth rates

            for K2N in K2_N:
                for K2I in K2_I:
                    for K2S in K2_S:                        # Drug kill rates
                        for ki in Deg:                      # Drug degradation rate
                            for kw in Ind:                  # Drug induction / input rate

                                for o in range(len(K3_NS)):
                                    # Interconversion rates between cell states
                                    K3NS = K3_NS[o]  # N -> S
                                    K3NI = K3_NI[o]  # N -> I
                                    K3IN = K3_IN[o]  # I -> N
                                    K3IS = K3_IS[o]  # I -> S
                                    K3SI = K3_SI[o]  # S -> I
                                    K3SN = K3_SN[o]  # S -> N

                                    row = {
                                        keys[0]: N,   keys[1]: I,   keys[2]: S,
                                        keys[3]: K1N, keys[4]: K1I, keys[5]: K1S,
                                        keys[6]: K2N, keys[7]: K2I, keys[8]: K2S,
                                        keys[9]:  K3NS, keys[10]: K3NI, keys[11]: K3IN,
                                        keys[12]: K3IS, keys[13]: K3SI, keys[14]: K3SN,
                                        keys[15]: ki, keys[16]: kw
                                    }
                                    data.append(row)

    # Save full parameter table to CSV for inspection
    csv_path = os.path.join(save_dir, 'InitialParamSpace.csv')
    df = pd.DataFrame(data)
    df.to_csv(csv_path, index=False)

    # Save each parameter combination as an individual pickle for parallel processing
    for i, row in enumerate(data, start=1):
        pickle_path = os.path.join(pklSave, f'Combo{i}.pkl')
        with open(pickle_path, 'wb') as f:
            pickle.dump(row, f)

    print("Files saved successfully!")


# =============================================================================
# BURN-IN MODEL  (no drug)
# =============================================================================
# Logistic growth with interconversion between three cell states (N, I, S).
# K  = carrying capacity (normalised to 1)
# K1 = intrinsic growth rate per state
# K3 = interconversion rate between states

def BurnInmodel(t, y, K, K1N, K1I, K1S,
                K3NS, K3NI, K3IN, K3IS, K3SI, K3SN):
    N, I, S = y

    dNdt = K1N * N * (1 - (N + I + S) / K) + K3SN * S + K3IN * I - K3NS * N - K3NI * N
    dIdt = K1I * I * (1 - (N + I + S) / K) + K3NI * N + K3SI * S - K3IS * I - K3IN * I
    dSdt = K1S * S * (1 - (N + I + S) / K) + K3NS * N + K3IS * I - K3SN * S - K3SI * S

    return [dNdt, dIdt, dSdt]


def termination_event(t, y, P_threshold):
    """Stop integration when total population reaches P_threshold."""
    return (y[0] + y[1] + y[2]) - P_threshold


# =============================================================================
# BURN-IN SIMULATION
# =============================================================================
# Grows the tumour from near-zero to a realistic starting population.
# For each base parameter set, 64 variants are run (all on/off combinations
# of the 6 K3 interconversion rates, encoded as a 6-bit binary number).

def BurnInSim(pkl_files, pklSave, result_dir, load_dir,
              time_steps, num_points, P_threshold,
              save_dir, BurnInFinal, data_base):

    start = pytime.time()
    final_values    = []
    time_series_data = []

    os.makedirs(save_dir,   exist_ok=True)
    os.makedirs(pklSave,    exist_ok=True)

    for pkl_file in pkl_files:

        # Load one base parameter combination
        file_path = os.path.join(load_dir, pkl_file)
        with open(file_path, 'rb') as f:
            base_data = pickle.load(f)

        K3_keys = ['K3NS', 'K3NI', 'K3IN', 'K3IS', 'K3SI', 'K3SN']

        # --- 64 binary variants: each bit controls whether a K3 rate is active ---
        for i in range(64):
            data   = base_data.copy()
            binary = format(i, '06b')   # e.g. '101010'

            # Switch each K3 on (original value) or off (0) per bit
            for j, key in enumerate(K3_keys):
                data[key] = base_data[key] if binary[j] == '1' else 0.0

            # Build a human-readable interconversion string, e.g. "<NI>|S|"
            interconversion = ""
            for j, key in enumerate(K3_keys):
                if data[key] == 0:
                    interconversion += "|"
                elif j % 2 == 0:
                    interconversion += "<"
                else:
                    interconversion += ">"
            interconversion = (interconversion[0:1] + 'N' + interconversion[1:2] +
                               interconversion[2:3] + 'I' + interconversion[3:4] +
                               interconversion[4:5] + 'S' + interconversion[5:6])

            initial_state = [data['N'], data['I'], data['S']]

            # Bind current parameters into the model function
            model_with_params = partial(
                BurnInmodel, K=1,
                K1N=data['K1N'], K1I=data['K1I'], K1S=data['K1S'],
                K3NS=data['K3NS'], K3NI=data['K3NI'],
                K3IN=data['K3IN'], K3IS=data['K3IS'],
                K3SI=data['K3SI'], K3SN=data['K3SN']
            )

            # Stop early if population exceeds P_threshold
            event_with_threshold = partial(termination_event, P_threshold=P_threshold)
            event_with_threshold.terminal = True

            sol = solve_ivp(
                fun=model_with_params,
                t_span=(0, time_steps),
                y0=initial_state,
                t_eval=np.linspace(0, time_steps, num_points),
                events=event_with_threshold,
                method='RK45'
            )

            Final_N = sol.y[0, -1]
            Final_I = sol.y[1, -1]
            Final_S = sol.y[2, -1]

            # --- Plot burn-in trajectory ---
            plt.figure(figsize=(10, 6))
            plt.plot(sol.t, sol.y[0, :], label="N Cell")
            plt.plot(sol.t, sol.y[1, :], label="I Cell")
            plt.plot(sol.t, sol.y[2, :], label="S Cell")
            plt.xlabel("Time")
            plt.ylabel("Population")
            plt.title("Time Evolution of N, I, S during Burn-In (NTG)")
            plt.legend()
            plt.grid(True)
            plot_file = os.path.join(result_dir, f"NTG_{pkl_file}.png")
            plt.savefig(plot_file)
            plt.close()
            print(f"Saved plot: {plot_file}")

            # Build RunID: base filename + variant index
            run_id = f"{pkl_file.replace('.pkl','')}_mod_{i}.pkl"

            row = {
                "RunID": run_id,
                "Interconversion": interconversion,
                "Final_N": Final_N, "Final_I": Final_I, "Final_S": Final_S,
                "kw": data['kw'], "ki": data['ki'],
                "K1N": data['K1N'], "K1I": data['K1I'], "K1S": data['K1S'],
                "K2N": data['K2N'], "K2I": data['K2I'], "K2S": data['K2S'],
                "K3NS": data['K3NS'], "K3NI": data['K3NI'], "K3IN": data['K3IN'],
                "K3IS": data['K3IS'], "K3SI": data['K3SI'], "K3SN": data['K3SN']
            }
            final_values.append(row)

            time_series_data.append({
                "RunID": run_id,
                "t_eval": sol.t.tolist(),
                "populations": list(zip(sol.y[0], sol.y[1], sol.y[2]))
            })

            # Save the final state of this variant as input for the treatment phase
            new_pkl_path = os.path.join(pklSave, run_id)
            with open(new_pkl_path, 'wb') as new_pkl:
                pickle.dump({
                    'Final_N': Final_N, 'Final_I': Final_I, 'Final_S': Final_S,
                    "K1N": data['K1N'], "K1I": data['K1I'], "K1S": data['K1S'],
                    "K2N": data['K2N'], "K2I": data['K2I'], "K2S": data['K2S'],
                    'K3NS': data['K3NS'], 'K3NI': data['K3NI'], 'K3IN': data['K3IN'],
                    'K3IS': data['K3IS'], 'K3SI': data['K3SI'], 'K3SN': data['K3SN'],
                    'kw': data['kw'], 'ki': data['ki'],
                }, new_pkl)

    # Save final summary to CSV
    pd.DataFrame(final_values).to_csv(os.path.join(save_dir, BurnInFinal), index=False)

    end = pytime.time()
    print(f"Burn-in complete. Took {end - start:.2f} seconds.")


# =============================================================================
# TREATMENT MODEL  (with drug)
# =============================================================================
# Extends the burn-in model with:
#   C    = drug concentration (4th state variable)
#   K2   = drug-induced kill rate per cell state
#   kw   = drug input/induction rate
#   ki   = drug degradation rate
# Drug kinetics: dC/dt = kw - ki*C  (simple first-order PK)

def Treamentmodel(t, y, K, K1N, K1I, K1S, K2N, K2I, K2S,
                  K3NS, K3NI, K3IN, K3IS, K3SI, K3SN,
                  kw, ki, C_param):
    N, I, S, C = y

    dNdt = (K1N * N * (1 - (N + I + S) / K)
            - K2N * N * C
            + K3SN * S + K3IN * I - K3NS * N - K3NI * N)

    dIdt = (K1I * I * (1 - (N + I + S) / K)
            - K2I * I * C
            + K3SI * S + K3NI * N - K3IS * I - K3IN * I)

    dSdt = (K1S * S * (1 - (N + I + S) / K)
            - K2S * S * C
            + K3NS * N + K3IS * I - K3SN * S - K3SI * S)

    dCdt = kw - ki * C     # First-order drug PK

    return [dNdt, dIdt, dSdt, dCdt]


def TreatmentSim(pklSave, pkl_files, load_dir, result_dir,
                 time_steps, num_points, P_threshold,
                 save_dir, TreatmentFinal, interval, data_base):

    start = pytime.time()
    final_values     = []
    time_series_data = []  # NOTE: was referenced before assignment in original code – initialised here

    if not pkl_files:
        print("No pickle files found in the directory.")
        return

    os.makedirs(save_dir,  exist_ok=True)
    os.makedirs(pklSave,   exist_ok=True)

    # Split the full time window into sub-intervals of length `interval`
    t_intervals = [[x, min(x + interval, time_steps)]
                   for x in range(0, time_steps, interval)]

    for pkl_file in pkl_files:

        file_path = os.path.join(load_dir, pkl_file)
        with open(file_path, 'rb') as f:
            data = pickle.load(f)
        print(data)

        C = 0.0  # Drug starts at zero at the beginning of treatment

        # Initial condition for treatment = final state from burn-in
        current_state = [data['Final_N'], data['Final_I'], data['Final_S'], C]

        t_total = []
        y_total = []

        K3_keys = ['K3NS', 'K3NI', 'K3IN', 'K3IS', 'K3SI', 'K3SN']
        K3s = {key: data[key] for key in K3_keys}

        # Rebuild interconversion label from K3 pattern
        interconversion = ""
        for j, key in enumerate(K3_keys):
            if K3s[key] == 0:
                interconversion += "|"
            elif j % 2 == 0:
                interconversion += "<"
            else:
                interconversion += ">"
        interconversion = (interconversion[0:1] + 'N' + interconversion[1:2] +
                           interconversion[2:3] + 'I' + interconversion[3:4] +
                           interconversion[4:5] + 'S' + interconversion[5:6])

        # --- Integrate piecewise over each sub-interval ---
        for timespan in t_intervals:
            t_evals = np.linspace(timespan[0], timespan[1], 101)

            # Flatten current_state to a 1-D float array (handles list/array inputs safely)
            y0 = np.asarray(
                [val[0] if isinstance(val, (list, np.ndarray)) else val
                 for val in current_state],
                dtype=np.float64
            ).flatten()

            sol = solve_ivp(
                fun=lambda t, y: Treamentmodel(
                    t, y, K=1,
                    K1N=data['K1N'], K1I=data['K1I'], K1S=data['K1S'],
                    K2N=data['K2N'], K2I=data['K2I'], K2S=data['K2S'],
                    K3NS=data['K3NS'], K3NI=data['K3NI'], K3IN=data['K3IN'],
                    K3IS=data['K3IS'], K3SI=data['K3SI'], K3SN=data['K3SN'],
                    kw=data['kw'], ki=data['ki'],
                    C_param=0
                ),
                t_span=timespan,
                y0=y0,
                method='RK45',
                t_eval=t_evals,
                vectorized=False,
                dense_output=True
            )

            t_total.append(sol.t)
            y_total.append(sol.y)

            # Pass the final state of this interval as the initial state for the next
            current_state = sol.y[:, -1]

        # Combine piecewise results into single arrays
        t_total = np.concatenate(t_total)
        y_total = np.concatenate(y_total, axis=1)

        time_series_data.append({
            "RunID": pkl_file,
            "t_eval": t_total.tolist(),
            "populations": list(zip(
                y_total[0].tolist(), y_total[1].tolist(),
                y_total[2].tolist(), y_total[3].tolist()
            ))
        })

        final_values.append({
            "RunID": pkl_file,
            'Interconversion': interconversion,
            'Initial_N': data['Final_N'], 'Initial_I': data['Final_I'], 'Initial_S': data['Final_S'],
            "Final_N": max(0, y_total[0, -1]),
            "Final_I": max(0, y_total[1, -1]),
            "Final_S": max(0, y_total[2, -1]),
            "Drug": data['kw'],
            "K1N": data['K1N'], "K1I": data['K1I'], "K1S": data['K1S'],
            "K2N": data['K2N'], "K2I": data['K2I'], "K2S": data['K2S'],
            "K3NS": data['K3NS'], "K3NI": data['K3NI'], "K3IN": data['K3IN'],
            "K3IS": data['K3IS'], "K3SI": data['K3SI'], "K3SN": data['K3SN']
        })

        # Save this run's final state for the relapse phase
        new_pkl_path = os.path.join(pklSave, pkl_file)
        with open(new_pkl_path, 'wb') as new_pkl:
            pickle.dump(final_values[-1], new_pkl)

        # --- Plot treatment trajectory ---
        plt.figure(figsize=(10, 6))
        plt.plot(t_total, y_total[0, :], label="N Cell")
        plt.plot(t_total, y_total[1, :], label="I Cell")
        plt.plot(t_total, y_total[2, :], label="S Cells")
        plt.plot(t_total, y_total[3, :], label="Drug Concentration")
        plt.xlabel("Time")
        plt.ylabel("Concentration / Population")
        plt.title("Time Evolution of N, I, S, and Drug during Treatment")
        plt.legend()
        plt.grid(True)
        plot_file = os.path.join(result_dir, f"Treatment_C_{pkl_file}.png")
        plt.savefig(plot_file)
        plt.close()
        print(f"Saved plot: {plot_file}")

    pd.DataFrame(final_values).to_csv(os.path.join(save_dir, TreatmentFinal), index=False)

    end = pytime.time()
    print(f"Treatment complete. Took {end - start:.2f} seconds.")


# =============================================================================
# RELAPSE MODEL  (no drug — mirrors burn-in ODE structure)
# =============================================================================
# After treatment ends, drug is removed and the remaining cells may regrow.
# *** FIX APPLIED HERE ***
# Populations are floored at 0 inside the ODE so the solver never uses
# negative values in its calculations.  Without this, the RK45 stepper can
# overshoot zero and then grow a biologically meaningless negative population.

def Relapsemodel(t, y, K1N, K1I, K1S,
                 K3NS, K3NI, K3IN, K3IS, K3SI, K3SN,
                 K):

    # Floor each population at zero – cells cannot be negative
    N = max(0.0, y[0])
    I = max(0.0, y[1])
    S = max(0.0, y[2])

    dNdt = (K1N * N * (1 - (N + I + S) / K)
            + K3SN * S + K3IN * I - K3NS * N - K3NI * N)

    dIdt = (K1I * I * (1 - (N + I + S) / K)
            + K3SI * S + K3NI * N - K3IS * I - K3IN * I)

    dSdt = (K1S * S * (1 - (N + I + S) / K)
            + K3NS * N + K3IS * I - K3SN * S - K3SI * S)

    return [dNdt, dIdt, dSdt]


def Relapse_termination(t, y, P_threshold):
    """
    Event function: triggers when total population (N+I+S) reaches P_threshold.
    Use floored values to be consistent with the ODE floor above.
    """
    N = max(0.0, y[0])
    I = max(0.0, y[1])
    S = max(0.0, y[2])
    return (N + I + S) - P_threshold


# =============================================================================
# RELAPSE SIMULATION
# =============================================================================
# Runs post-treatment regrowth starting from the final treatment state.
# Stops early if total population reaches P_threshold (relapse criterion).

def RelapseSim(pkl_files, load_dir, time_steps, num_points,
               P_threshold, save_dir, RelapseFinal, data_base):

    final_values = []
    start = pytime.time()   # BUG FIX: original code called pytime.time (no call brackets)

    if not pkl_files:
        print("No pickle files found in the directory.")
        return

    os.makedirs(save_dir, exist_ok=True)

    for pkl_file in pkl_files:

        file_path = os.path.join(load_dir, pkl_file)

        if not os.path.exists(file_path):
            print(f"File not found: {file_path}")
            continue

        with open(file_path, 'rb') as f:
            data = pickle.load(f)

        # Starting populations = final state after treatment
        initial_state = [data['Final_N'], data['Final_I'], data['Final_S']]

        K1N, K1I, K1S   = data['K1N'], data['K1I'], data['K1S']
        K2N, K2I, K2S   = data['K2N'], data['K2I'], data['K2S']
        Drug             = data['Drug']
        interconversion  = data['Interconversion']

        K3_keys = ['K3NS', 'K3NI', 'K3IN', 'K3IS', 'K3SI', 'K3SN']
        K3s = {key: data[key] for key in K3_keys}

        # Bind parameters into the model
        model_with_params = partial(
            Relapsemodel, K=1,
            K1N=K1N, K1I=K1I, K1S=K1S,
            K3NS=K3s['K3NS'], K3NI=K3s['K3NI'], K3IN=K3s['K3IN'],
            K3IS=K3s['K3IS'], K3SI=K3s['K3SI'], K3SN=K3s['K3SN']
        )

        # Stop if total population reaches P_threshold (relapse)
        event_with_threshold = partial(Relapse_termination, P_threshold=P_threshold)
        event_with_threshold.terminal  = True
        event_with_threshold.direction = 0   # trigger on both rising and falling edges

        t_span = (0, time_steps)
        t_eval = np.linspace(0, time_steps, num_points)

        sol = integrate.solve_ivp(
            fun=model_with_params,
            t_span=t_span,
            y0=initial_state,
            method='RK45',
            t_eval=t_eval,
            events=event_with_threshold,
            vectorized=False
        )

        # Floor final values – belt-and-braces guard on the output
        Final_N = max(0.0, sol.y[0, -1])
        Final_I = max(0.0, sol.y[1, -1])
        Final_S = max(0.0, sol.y[2, -1])

        row = {
            "RunID":          pkl_file,
            'Initial_N':      data['Final_N'],
            'Initial_I':      data['Final_I'],
            'Initial_S':      data['Final_S'],
            'Interconversion': interconversion,
            'Drug':           Drug,
            'K1N': K1N,  'K1I': K1I,  'K1S': K1S,
            'K3NS': K3s['K3NS'], 'K3NI': K3s['K3NI'], 'K3IN': K3s['K3IN'],
            'K3IS': K3s['K3IS'], 'K3SI': K3s['K3SI'], 'K3SN': K3s['K3SN'],
            'Final_N': Final_N, 'Final_I': Final_I, 'Final_S': Final_S
        }
        final_values.append(row)

    # Save relapse summary to CSV
    csv_path = os.path.join(save_dir, RelapseFinal)
    pd.DataFrame(final_values).to_csv(csv_path, index=False)

    end = pytime.time()
    print(f"Relapse complete. Took {end - start:.2f} seconds.")


# =============================================================================
# ENTRY POINTS – change paths below to match your directory structure
# =============================================================================

# --- Step 1: Generate parameter space ---
save_dir = r'C:\users\892215\OneDrive - hull.ac.uk\Desktop\Year3\Experiment_3\InSilico\S_Drug\Data'
pklSave  = r'C:\users\892215\OneDrive - hull.ac.uk\Desktop\Year3\Experiment_3\InSilico\S_Drug\BurnInIC'
os.makedirs(save_dir, exist_ok=True)

# Initial cell populations (near-zero: tumour seeding)
n = [0.000001]
m = [0.000001]
s = [0.000001]

# Intrinsic growth rates (N fastest, S slowest – resistant hierarchy)
K1_N = [0.75]
K1_I = [0.5]
K1_S = [0.25]

# Drug kill rates (K2S=1 always; K2I and K2N can be 0 or 0.5)
K2_N = [0, 0.5]
K2_I = [0, 0.5]
K2_S = [1]

# Interconversion rate values (applied combinatorially via binary switching in BurnInSim)
K3_NS = [0.1, 0.5, 0.9]
K3_NI = [0.1, 0.5, 0.9]
K3_IN = [0.1, 0.5, 0.9]
K3_IS = [0.1, 0.5, 0.9]
K3_SI = [0.1, 0.5, 0.9]
K3_SN = [0.1, 0.5, 0.9]

# Drug PK parameters: Concentration at steady state = Ind / Deg
Ind = [0.05, 0.25, 0.45]   # Drug input rates  -> steady-state C = {0.1, 0.5, 0.9}
Deg = [0.5]                 # Drug degradation rate

ParamSpace(pklSave, save_dir, n, m, s,
           K1_N, K1_I, K1_S,
           K2_N, K2_I, K2_S,
           K3_NS, K3_NI, K3_IN, K3_IS, K3_SI, K3_SN,
           Deg, Ind)

# --- Step 2: Burn-in simulation ---
load_dir   = r'C:\users\892215\OneDrive - hull.ac.uk\Desktop\Year3\Experiment_3\InSilico\S_Drug\BurnInIC'
save_dir   = r'C:\users\892215\OneDrive - hull.ac.uk\Desktop\Year3\Experiment_3\InSilico\S_Drug\Data'
pklSave    = r'C:\users\892215\OneDrive - hull.ac.uk\Desktop\Year3\Experiment_3\InSilico\S_Drug\TreatIC'
result_dir = r'C:\users\892215\OneDrive - hull.ac.uk\Desktop\Year3\Experiment_3\InSilico\S_Drug\Results'

data_base   = 'Time_Series.db'
BurnInFinal = 'BurnInFinal.csv'
os.makedirs(load_dir, exist_ok=True)

pkl_files        = [f for f in os.listdir(load_dir) if f.endswith('.pkl')]
final_values     = []
time_series_data = []

P_threshold = 0.4
time_steps  = 100
num_points  = 101

BurnInSim(pkl_files, pklSave, result_dir, load_dir,
          time_steps, num_points, P_threshold,
          save_dir, BurnInFinal, data_base)

# --- Step 3: Treatment simulation ---
load_dir   = r'C:\users\892215\OneDrive - hull.ac.uk\Desktop\Year3\Experiment_3\InSilico\S_Drug\TreatIC'
save_dir   = r'C:\users\892215\OneDrive - hull.ac.uk\Desktop\Year3\Experiment_3\InSilico\S_Drug\Data'
pklSave    = r'C:\users\892215\OneDrive - hull.ac.uk\Desktop\Year3\Experiment_3\InSilico\S_Drug\RelapsIC'
result_dir = r'C:\users\892215\OneDrive - hull.ac.uk\Desktop\Year3\Experiment_3\InSilico\S_Drug\Results'

data_base       = 'Time_Series.db'
TreatmentFinal  = 'TreatmentFinal.csv'
os.makedirs(load_dir, exist_ok=True)

pkl_files        = [f for f in os.listdir(load_dir) if f.endswith('.pkl')]
final_values     = []
time_series_data = []

interval    = 1
P_threshold = 0
time_steps  = 100
num_points  = 101

TreatmentSim(pklSave, pkl_files, load_dir, result_dir,
             time_steps, num_points, P_threshold,
             save_dir, TreatmentFinal, interval, data_base)

# --- Step 4: Relapse simulation ---
load_dir = r'C:\users\892215\OneDrive - hull.ac.uk\Desktop\Year3\Experiment_3\InSilico\S_Drug\RelapsIC'
save_dir = r'C:\users\892215\OneDrive - hull.ac.uk\Desktop\Year3\Experiment_3\InSilico\S_Drug\Data'

data_base    = 'Time_Series.db'
RelapseFinal = 'RelapseFinal.csv'
os.makedirs(load_dir, exist_ok=True)

pkl_files        = [f for f in os.listdir(load_dir) if f.endswith('.pkl')]
final_values     = []
time_series_data = []

P_threshold = 0.9
time_steps  = 100
num_points  = 101

RelapseSim(pkl_files, load_dir, time_steps, num_points,
           P_threshold, save_dir, RelapseFinal, data_base)
