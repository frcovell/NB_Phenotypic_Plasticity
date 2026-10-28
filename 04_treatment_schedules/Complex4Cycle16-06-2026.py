# Schedule D (four-phase) used for Paper 1, Fig. 11.
# Reconstructed from Complex16-06-2026.py (3-cycle version): only two lines changed -
#   save folder "3 cycle" -> "4 cycle" and cycle_length = 3 -> 4.
import numpy as np
import pandas as pd
from scipy.stats import qmc
from scipy.integrate import solve_ivp 
from functools import partial
import os
import matplotlib.pyplot as plt
import time
import sqlite3





multi_complex_save_dir = r'C:\Users\892215\OneDrive - hull.ac.uk\Desktop\Year3\Experiment_3\Validity\Plots\4 cycle\Multi ComplexTS'

#os.makedirs(plot_dir, exist_ok=True)
# -------------------------- INITIAL CONDITIONS --------------------------

#set up latin hypercube, 100 samples of 3 values
sampler = qmc.LatinHypercube(d=3, seed=42)
sample = sampler.random(n=10000)

# Normalise each row to sum to 1 then scale to 0.4
sample_normalised = sample / sample.sum(axis=1, keepdims=True) * 0.4
valid = pd.DataFrame(sample_normalised, columns=['Initial N', 'Initial I', 'Initial S'])
#print(valid)
# True dominant cases - cell must be greater than BOTH others
n_dom_candidates = valid[(valid['Initial N'] > valid['Initial I']) & 
                          (valid['Initial N'] > valid['Initial S'])]
i_dom_candidates = valid[(valid['Initial I'] > valid['Initial N']) & 
                          (valid['Initial I'] > valid['Initial S'])]
s_dom_candidates = valid[(valid['Initial S'] > valid['Initial N']) & 
                          (valid['Initial S'] > valid['Initial I'])]
rounded = valid.round(3)
no_dom_candidates = rounded[(rounded['Initial S'] == rounded['Initial N']) | 
                          (rounded['Initial S'] == rounded['Initial I']) |
                          (rounded['Initial N'] == rounded['Initial I'])]

# Get top 5 most dominant cases for each
n_dom_idx = n_dom_candidates.nlargest(50, 'Initial N').index
i_dom_idx = i_dom_candidates.nlargest(50, 'Initial I').index
s_dom_idx = s_dom_candidates.nlargest(50, 'Initial S').index
no_dom_idx = no_dom_candidates.head(49)#nlargest(50).index
#print(type(no_dom_idx))
no_dom = pd.concat([
    no_dom_idx,
    pd.DataFrame([{'Initial N': 0.4/3, 'Initial I': 0.4/3, 'Initial S': 0.4/3}])
], ignore_index=True)
print(no_dom)

InitialPop = pd.concat([
    no_dom.assign(Dominant='No Dominant'),
    valid.loc[n_dom_idx].assign(Dominant='N Dominant'),
    valid.loc[i_dom_idx].assign(Dominant='I Dominant'),
    valid.loc[s_dom_idx].assign(Dominant='S Dominant'),
]).reset_index(drop=True)

# Round then calculate total
InitialPop = InitialPop.round(6)
InitialPop['Total'] = InitialPop[['Initial N', 'Initial I', 'Initial S']].sum(axis=1)

print(InitialPop)
print(f"Total initial conditions: {len(InitialPop)}")
# -------------------------- PARAMETERS --------------------------
# Growth rates
K1_N = 0.75
K1_I = 0.5
K1_S = 0.25

# Death rates Phase 2
Ph2K2_N = 1  # target cell
Ph2K2_I=[0,0.5] # non target
Ph2K2_S=[0,0.5]  # non target

# Death rates Phase 4
Ph4K2_S = 1  # target cell
Ph4K2_I=[0,0.5] # non target
Ph4K2_N=[0,0.5]  # non target


# Drug
kw = 0.4
ki = 0.5

# Phase 1 and 2 patterns
k3_combinations = [
    (0, 0, 1, 1, 0, 1),  # P3 |N| <I> |S>
    (0, 0, 1, 1, 0, 0),  # P4 |N| <I> |S|
    (0, 0, 1, 1, 1, 1)   # P5 |N| <I> <S>
]
pattern_names = {
    (0, 0, 1, 1, 0, 1): 'P3',
    (0, 0, 1, 1, 0, 0): 'P4',
    (0, 0, 1, 1, 1, 1): 'P5'
}

P1 = (1,1,0,1,0,0) #<N> |I> |S|

# Interconversion rates
K3_rates = [0.5]

#number of cycles
cycle_length = 4
# -------------------------- Model Function --------------------------
# Define the model function
def Treamentmodel(t, y, K, K1N, K1I, K1S, K2N, K2I, K2S,
                  K3NS, K3NI, K3IN, K3IS, K3SI, K3SN,
                  kw, ki, C_param):  # Avoid overwriting C
    N, I, S, C = y  # Ensure correct unpacking
    
    # Define ODEs
    dNdt = K1N * N * (1 - (N + I + S) / K) - K2N * N * C + K3SN * S + K3IN * I - K3NS * N - K3NI * N
    dIdt = K1I * I * (1 - (N + I + S) / K) - K2I * I * C + K3SI * S + K3NI * N - K3IS * I - K3IN * I
    dSdt = K1S * S * (1 - (N + I + S) / K) - K2S * S * C + K3NS * N + K3IS * I - K3SN * S - K3SI * S
    dCdt = kw - ki * C
    
    return [dNdt, dIdt, dSdt, dCdt]

# -------------------------- THRESHOLD HELPERS --------------------------
# setting a threshold of 1e-6 and capture when the total population drops below 
def find_threshold_time(sol, threshold=1e-3):
    """Find first timestep where total population drops to or below threshold"""
    #save total population at every time step
    total_pop = sol.y[0, :] + sol.y[1, :] + sol.y[2, :]
    # Check which time step total population chrosses threshold 
    below_threshold = np.where(total_pop <= threshold)[0]
    # If a time step is recorded save it
    if len(below_threshold) > 0:
        return float(sol.t[below_threshold[0]])
    elif total_pop[-1] <= threshold:
        return float(sol.t[-1])  # return final timestep if end value is below threshold
    else:
        return None
        
# setting threshold for relapse to check if populations rises above 1e-6
def find_recovery_time(sol, threshold=1e-3):
    """Find first timestep where total population rises above threshold"""
    total_pop = sol.y[0, :] + sol.y[1, :] + sol.y[2, :]
    above_threshold = np.where(total_pop > threshold)[0]
    if len(above_threshold) > 0:
        return float(sol.t[above_threshold[0]])
    elif total_pop[-1] > threshold:
        return float(sol.t[-1])  # return final timestep if end value is above threshold
    else:
        return None
# -------------------------- connect to database --------------------------

# Connect BEFORE the loop
os.makedirs(multi_complex_save_dir, exist_ok=True)
conn = sqlite3.connect(os.path.join(multi_complex_save_dir, 'timeseries.db'))
cur = conn.cursor()
cur.execute('DROP TABLE IF EXISTS Complex_TimeSeries')
cur.execute('''CREATE TABLE Complex_TimeSeries (
    dominant    TEXT,
    pattern_id  TEXT,
    k3_rate     REAL,
    ph2k2i      REAL,
    ph2k2s      REAL,
    ph4k2i      REAL,
    ph4k2n      REAL,
    cycle       INTEGER,
    phase       INTEGER,
    time        REAL,
    N           REAL,
    I           REAL,
    S           REAL
)''')
conn.commit()
    # -------------------------- SIMULATION --------------------------
results = []
dominant_counts = {}

for i in range(len(InitialPop)):
    N0 = InitialPop.iloc[i].loc['Initial N']
    I0 = InitialPop.iloc[i].loc['Initial I']
    S0 = InitialPop.iloc[i].loc['Initial S']
    dominant = InitialPop.iloc[i].loc['Dominant']

    # Track count per dominant type
    if dominant not in dominant_counts:
        dominant_counts[dominant] = 1
    else:
        dominant_counts[dominant] += 1
    dom_idx = dominant_counts[dominant]

    for Ph2K2I in Ph2K2_I:
        for Ph2K2S in Ph2K2_S:
            for Ph4K2I in Ph4K2_I:
                for Ph4K2N in Ph4K2_N:
                    for pattern in k3_combinations:
                        pattern_name = pattern_names[pattern]
                        K3NS, K3NI, K3IN, K3IS, K3SI, K3SN = pattern
                        P1K3NS, P1K3NI, P1K3IN, P1K3IS, P1K3SI, P1K3SN = P1

                        K3_keys = ['K3NS', 'K3NI', 'K3IN', 'K3IS', 'K3SI', 'K3SN']
                        K3_vals = [K3NS, K3NI, K3IN, K3IS, K3SI, K3SN]

                        interconversion = ""
                        for j, val in enumerate(K3_vals):
                            if val == 0:
                                interconversion += "|"
                            elif j % 2 == 0:
                                interconversion += "<"
                            else:
                                interconversion += ">"
                        interconversion = (interconversion[0:1] + 'N' + interconversion[1:2] +
                                          interconversion[2:3] + 'I' + interconversion[3:4] +
                                          interconversion[4:5] + 'S' + interconversion[5:6])

                        for K3_rate in K3_rates:
                            K3NS_val = float(K3NS * K3_rate)
                            K3NI_val = float(K3NI * K3_rate)
                            K3IN_val = float(K3IN * K3_rate)
                            K3IS_val = float(K3IS * K3_rate)
                            K3SI_val = float(K3SI * K3_rate)
                            K3SN_val = float(K3SN * K3_rate)

                            P1K3NS_val = float(P1K3NS * K3_rate)
                            P1K3NI_val = float(P1K3NI * K3_rate)
                            P1K3IN_val = float(P1K3IN * K3_rate)
                            P1K3IS_val = float(P1K3IS * K3_rate)
                            P1K3SI_val = float(P1K3SI * K3_rate)
                            P1K3SN_val = float(P1K3SN * K3_rate)

                            cycle_solutions = []
                            current_y0 = [N0, I0, S0, 0]
                            # Reset before cycle loop
                            prev_phase4_threshold_time = None

                            for q in range(cycle_length):
                                cycle_start = time.time()
                               # Phase 1 - P3/P4/P5 pattern induction, no drug
                                model_phase1 = partial(Treamentmodel, K=1,
                                    K1N=K1_N, K1I=K1_I, K1S=K1_S,
                                    K2N=0, K2I=0, K2S=0,
                                    K3NS=K3NS_val, K3NI=K3NI_val, K3IN=K3IN_val,
                                    K3IS=K3IS_val, K3SI=K3SI_val, K3SN=K3SN_val,
                                    kw=0, ki=ki, C_param=0)

                                sol_phase1 = solve_ivp(
                                    fun=model_phase1,
                                    t_span=[0, 30],
                                    y0=current_y0,
                                    method='RK45',
                                    t_eval=np.linspace(0, 30, 31)
                                )

                                # Phase 2 - N inhibitor
                                model_phase2 = partial(Treamentmodel, K=1,
                                    K1N=K1_N, K1I=K1_I, K1S=K1_S,
                                    K2N=Ph2K2_N, K2I=float(Ph2K2I), K2S=float(Ph2K2S),
                                    K3NS=K3NS_val, K3NI=K3NI_val, K3IN=K3IN_val,
                                    K3IS=K3IS_val, K3SI=K3SI_val, K3SN=K3SN_val,
                                    kw=kw, ki=ki, C_param=0)

                                sol_phase2 = solve_ivp(
                                    fun=model_phase2,
                                    t_span=[0, 100],
                                    y0=sol_phase1.y[:, -1],
                                    method='RK45',
                                    t_eval=np.linspace(0, 100, 101)
                                )

                                # Phase 3 - P1 pattern induction, no drug
                                model_phase3 = partial(Treamentmodel, K=1,
                                    K1N=K1_N, K1I=K1_I, K1S=K1_S,
                                    K2N=0, K2I=0, K2S=0,
                                    K3NS=P1K3NS_val, K3NI=P1K3NI_val, K3IN=P1K3IN_val,
                                    K3IS=P1K3IS_val, K3SI=P1K3SI_val, K3SN=P1K3SN_val,
                                    kw=0, ki=ki, C_param=0)

                                sol_phase3 = solve_ivp(
                                    fun=model_phase3,
                                    t_span=[0, 30],
                                    y0=sol_phase2.y[:, -1],
                                    method='RK45',
                                    t_eval=np.linspace(0, 30, 31)
                                )

                                # Phase 4 - S inhibitor
                                model_phase4 = partial(Treamentmodel, K=1,
                                    K1N=K1_N, K1I=K1_I, K1S=K1_S,
                                    K2N=float(Ph4K2N), K2I=float(Ph4K2I), K2S=Ph4K2_S,
                                    K3NS=P1K3NS_val, K3NI=P1K3NI_val, K3IN=P1K3IN_val,
                                    K3IS=P1K3IS_val, K3SI=P1K3SI_val, K3SN=P1K3SN_val,
                                    kw=kw, ki=ki, C_param=0)

                                sol_phase4 = solve_ivp(
                                    fun=model_phase4,
                                    t_span=[0, 100],
                                    y0=sol_phase3.y[:, -1],
                                    method='RK45',
                                    t_eval=np.linspace(0, 100, 101)
                                )

                                cycle_solutions.append({
                                    'cycle': q + 1,
                                    'sol_phase1': sol_phase1,
                                    'sol_phase2': sol_phase2,
                                    'sol_phase3': sol_phase3,
                                    'sol_phase4': sol_phase4
                                })

                                current_y0 = np.clip(sol_phase4.y[:, -1], 0, None)

                                # Calculate final totals
                                phase1_total = (max(0, sol_phase1.y[0, -1]) +
                                               max(0, sol_phase1.y[1, -1]) +
                                               max(0, sol_phase1.y[2, -1]))
                                phase2_total = (max(0, sol_phase2.y[0, -1]) +
                                               max(0, sol_phase2.y[1, -1]) +
                                               max(0, sol_phase2.y[2, -1]))
                                phase3_total = (max(0, sol_phase3.y[0, -1]) +
                                               max(0, sol_phase3.y[1, -1]) +
                                               max(0, sol_phase3.y[2, -1]))
                                phase4_total = (max(0, sol_phase4.y[0, -1]) +
                                               max(0, sol_phase4.y[1, -1]) +
                                               max(0, sol_phase4.y[2, -1]))

                                # Find threshold and recovery times
                                phase1_recovery_time = find_recovery_time(sol_phase1)
                                phase2_threshold_time = find_threshold_time(sol_phase2, threshold= 1e-3)
                                phase3_recovery_time = find_recovery_time(sol_phase3)
                                phase4_threshold_time = find_threshold_time(sol_phase4, threshold=1e-3)
                                phase4_recovery_time = find_recovery_time(sol_phase4)

                                results.append({
                                    'Dominant': dominant,
                                    'Initial N': N0,
                                    'Initial I': I0,
                                    'Initial S': S0,
                                    'Ph2K2I': Ph2K2I,
                                    'Ph2K2S': Ph2K2S,
                                    'Ph4K2I': Ph4K2I,
                                    'Ph4K2N': Ph4K2N,
                                    'Pattern ID': pattern_name,
                                    'Pattern': interconversion,
                                    'K3_rate': K3_rate,
                                    'Cycle': q + 1,
                                    'Phase1_Final_N': max(0, sol_phase1.y[0, -1]),
                                    'Phase1_Final_I': max(0, sol_phase1.y[1, -1]),
                                    'Phase1_Final_S': max(0, sol_phase1.y[2, -1]),
                                    'Phase1_Total': phase1_total,
                                     # Use PREVIOUS cycle's phase4_threshold_time
                                    # Phase1_Relapse in results.append:
                                    'Phase1_Relapse': (q > 0) and prev_phase4_threshold_time is not None and phase1_total > 1e-3,
                                    'Phase1_Recovery_Time': phase1_recovery_time if (q > 0 and prev_phase4_threshold_time is not None) else None,

                                    'Phase2_Final_N': max(0, sol_phase2.y[0, -1]),
                                    'Phase2_Final_I': max(0, sol_phase2.y[1, -1]),
                                    'Phase2_Final_S': max(0, sol_phase2.y[2, -1]),
                                    'Phase2_Total': phase2_total,
                                    'Phase2_Reached_Zero': phase2_threshold_time is not None and phase2_total <= 1e-3,
                                    'Phase2_Zero_Time': phase2_threshold_time,
                                    'Phase3_Final_N': max(0, sol_phase3.y[0, -1]),
                                    'Phase3_Final_I': max(0, sol_phase3.y[1, -1]),
                                    'Phase3_Final_S': max(0, sol_phase3.y[2, -1]),
                                    'Phase3_Total': phase3_total,
                                    # check is remission was achived in phase 2 then cross threshold.
                                    'Phase3_Relapse': phase2_threshold_time is not None and phase3_total > 1e-3,
                                    # save time to cross threshold is phase 2 went into remission
                                    'Phase3_Recovery_Time': phase3_recovery_time if phase2_threshold_time is not None else None,
                                    'Phase4_Final_N': max(0, sol_phase4.y[0, -1]),
                                    'Phase4_Final_I': max(0, sol_phase4.y[1, -1]),
                                    'Phase4_Final_S': max(0, sol_phase4.y[2, -1]),
                                    'Phase4_Total': phase4_total,
                                    'Phase4_Reached_Zero': phase4_threshold_time is not None and phase4_total <= 1e-3,
                                    'Phase4_Zero_Time': phase4_threshold_time,
                                    'Phase4_Relapse': phase3_total > 1e-3 and phase4_threshold_time is not None and phase4_total <= 1e-3,
                                    'Phase4_Recovery_Time': phase4_recovery_time if phase4_threshold_time is not None else None,
                                })
                                
                                # Calculate phase4_threshold_time AFTER append so it's available for NEXT cycle's Phase1_Relapse
                                #phase4_threshold_time = find_threshold_time(sol_phase4, threshold= 1e-3)
                                # Replace the redundant line at the bottom of the loop:
                                prev_phase4_threshold_time = phase4_threshold_time

                                cycle_time = time.time() - cycle_start
                                print(f"Cycle {q+1} done | {dominant} | {pattern_name} | K3={K3_rate} | {cycle_time:.1f}s")
                                for phase_num, sol in [(1, sol_phase1), (2, sol_phase2), (3, sol_phase3), (4, sol_phase4)]:
                                    rows = [
                                        (dominant, pattern_name, K3_rate, Ph2K2I,Ph2K2S,Ph4K2I,Ph4K2N,
                                        q + 1, phase_num,
                                        float(sol.t[t_idx]),
                                        float(max(0, sol.y[0, t_idx])),
                                        float(max(0, sol.y[1, t_idx])),
                                        float(max(0, sol.y[2, t_idx])))
                                        for t_idx in range(len(sol.t))
                                    ]
                                    cur.executemany('''INSERT INTO Complex_TimeSeries
                                        (dominant, pattern_id, k3_rate, ph2k2i, ph2k2s, ph4k2i, ph4k2n,
                                        cycle, phase, time, N, I, S)
                                        VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)''', rows)

                                conn.commit()  
                            # Plot all cycles after all cycles complete
                            #fig, axes = plt.subplots(cycle_length, 4, figsize=(20, 5 * cycle_length))
#                            fig, axes = plt.subplots(cycle_length, 4, figsize=(16, 4 * cycle_length))
 #                           for q, cyc in enumerate(cycle_solutions):
  #                              t_phase1 = cyc['sol_phase1'].t
   #                             if t_phase1[-1] < 30:
    #                                print(f"Warning: Phase 1 solver stopped early at t={t_phase1[-1]:.2f} for cycle {cyc['cycle']}, {dominant}, {pattern_name}, K3={K3_rate}")
     #                           axes[q, 0].plot(t_phase1, cyc['sol_phase1'].y[0, :], label="N Cell", color='blue')
      #                          axes[q, 0].plot(t_phase1, cyc['sol_phase1'].y[1, :], label="I Cell", color='green')
       #                         axes[q, 0].plot(t_phase1, cyc['sol_phase1'].y[2, :], label="S Cell", color='red')
        #                        axes[q, 0].set_title(f'Cycle {cyc["cycle"]} Phase 1 - {pattern_name} Induction')
         #                       axes[q, 0].set_ylim(bottom=0)
          #                      axes[q, 0].set_xlabel('Time')
           #                     axes[q, 0].set_ylabel('Population')
            #                    axes[q, 0].legend()
             #                   axes[q, 0].grid(True)

#                                t_phase2 = cyc['sol_phase2'].t
 #                               if t_phase2[-1] < 100:
  #                                  print(f"Warning: Phase 2 solver stopped early at t={t_phase2[-1]:.2f} for cycle {cyc['cycle']}, {dominant}, {pattern_name}, K3={K3_rate}")
   #                             axes[q, 1].plot(t_phase2, cyc['sol_phase2'].y[0, :], label="N Cell", color='blue')
    #                            axes[q, 1].plot(t_phase2, cyc['sol_phase2'].y[1, :], label="I Cell", color='green')
     #                           axes[q, 1].plot(t_phase2, cyc['sol_phase2'].y[2, :], label="S Cell", color='red')
      #                          axes[q, 1].plot(t_phase2, cyc['sol_phase2'].y[3, :], label="Drug", color='orange', linestyle='--')
       #                         axes[q, 1].set_title(f'Cycle {cyc["cycle"]} Phase 2 - N Inhibitor\nK2I={Ph2K2I}, K2S={Ph2K2S}')
        ##                        axes[q, 1].set_ylim(bottom=0)
           #                     axes[q, 1].set_xlabel('Time')
          #                      axes[q, 1].set_ylabel('Population')
            #                    axes[q, 1].legend()
             #                   axes[q, 1].grid(True)

#                                t_phase3 = cyc['sol_phase3'].t
 #                               if t_phase3[-1] < 30:
  #                                  print(f"Warning: Phase 3 solver stopped early at t={t_phase3[-1]:.2f} for cycle {cyc['cycle']}, {dominant}, {pattern_name}, K3={K3_rate}")
   #                             axes[q, 2].plot(t_phase3, cyc['sol_phase3'].y[0, :], label="N Cell", color='blue')
    #                            axes[q, 2].plot(t_phase3, cyc['sol_phase3'].y[1, :], label="I Cell", color='green')
     #                           axes[q, 2].plot(t_phase3, cyc['sol_phase3'].y[2, :], label="S Cell", color='red')
         #                       axes[q, 2].set_title(f'Cycle {cyc["cycle"]} Phase 3 - P1 Induction')
      #                          axes[q, 2].set_ylim(bottom=0)
       #                         axes[q, 2].set_xlabel('Time')
        #                        axes[q, 2].set_ylabel('Population')
          #                      axes[q, 2].legend()
           #                     axes[q, 2].grid(True)

            #                    t_phase4 = cyc['sol_phase4'].t
             #                   if t_phase4[-1] < 100:
              #                      print(f"Warning: Phase 4 solver stopped early at t={t_phase4[-1]:.2f} for cycle {cyc['cycle']}, {dominant}, {pattern_name}, K3={K3_rate}")
               #                 axes[q, 3].plot(t_phase4, cyc['sol_phase4'].y[0, :], label="N Cell", color='blue')
                #                axes[q, 3].plot(t_phase4, cyc['sol_phase4'].y[1, :], label="I Cell", color='green')
                 #               axes[q, 3].plot(t_phase4, cyc['sol_phase4'].y[2, :], label="S Cell", color='red')
                  #              axes[q, 3].plot(t_phase4, cyc['sol_phase4'].y[3, :], label="Drug", color='orange', linestyle='--')
                   #             axes[q, 3].set_title(f'Cycle {cyc["cycle"]} Phase 4 - S Inhibitor\nK2I={Ph4K2I}, K2N={Ph4K2N}')
                    #            axes[q, 3].set_ylim(bottom=0)
                     #           axes[q, 3].set_xlabel('Time')
                      #          axes[q, 3].set_ylabel('Population')
                       #         axes[q, 3].legend()
                        #        axes[q, 3].grid(True)
#
 #                           plt.suptitle(f'Dominant: {dominant} | Pattern: {pattern_name} | K3={K3_rate} | Ph2K2I={Ph2K2I} | Ph2K2S={Ph2K2S} | Ph4K2I={Ph4K2I} | Ph4K2N={Ph4K2N}', fontsize=12)
  #                          plt.tight_layout()

                           # plot_file = os.path.join(multi_complex_save_dir, f'{dominant}_{pattern_name}_Ph2K2I{Ph2K2I}_Ph4K2I{Ph4K2I}_K3{K3_rate}.png')
   #                         plot_file = os.path.join(multi_complex_save_dir, f'{dominant}{dom_idx}_{pattern_name}_Ph2K2I{Ph2K2I}_Ph4K2I{Ph4K2I}_K3{K3_rate}.png')
#
 #                           plt.savefig(plot_file, dpi=72, bbox_inches='tight')
  #                          plt.close()
cur.execute('''CREATE INDEX idx_dom_pat 
               ON Complex_TimeSeries (dominant, pattern_id, cycle, phase)''')
conn.commit()
conn.close()
results_df = pd.DataFrame(results)
results_df.to_csv(os.path.join(multi_complex_save_dir, 'ComplexResults.csv'), index=False)
print(results_df)
print(f"Total runs: {len(results_df)}")