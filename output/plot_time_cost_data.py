#%%
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# ============================================
# 1. Load monomial-level CSV
# ============================================

#root = "beam_o_11_eps_11_Fmode_1_2026-05-19T11_43_07"
root = "beam_o_9_eps_9_Fmode_1_2026-05-20T01_11_38"
mono_file = root + "/benchmark_per_monomial.csv"
df_mono = pd.read_csv(mono_file, na_values='NaN')
df_mono['monomial_total_time_s'] = df_mono['monomial_total_time_s'].fillna(0)

# Global monomial index (1‑based row number)
global_index = np.arange(1, len(df_mono) + 1)
cumulative_mono = df_mono['monomial_total_time_s'].cumsum().values

# ============================================
# 2. Load order-level CSV
# ============================================
order_file = root + "/benchmark_per_order.csv"
df_order = pd.read_csv(order_file)

# Total time per order (sum of all three contributions)
df_order['total_time_order'] = (
    df_order['fillrhsG_time_s'] +
    df_order['fillrhsH_time_s'] +
    df_order['fillWf_time_s']
)
cumulative_order_totals = df_order['total_time_order'].cumsum().values

# Compute global start index of each order (first monomial index)
df_order['cumulative_monomials'] = df_order['n_monomials'].cumsum().values
start_indices = [1] + (df_order['cumulative_monomials'][:-1] + 1).tolist()

# ============================================
# 3. Build the step function for cumulative order cost
#    For each monomial, assign the cumulative order total of its order
# ============================================
cum_order_step = np.zeros(len(df_mono))
for i, start in enumerate(start_indices):
    end = df_order['cumulative_monomials'][i]          # last monomial of this order
    cum_val = cumulative_order_totals[i]
    cum_order_step[start-1:end] = cum_val             # fill from start to end (0‑based)

# Total cumulative = monomial cumulative + order step
total_cumulative = cumulative_mono + cum_order_step

# ============================================
# 3.5. Polynomial interpolation at order endpoints
# ============================================
#%%
import numpy as np
from numpy.polynomial import Polynomial
end_indices = df_order['cumulative_monomials'].tolist()  # last monomial of each order
# Use the end indices (global index of last monomial of each order)
x_endpoints = np.array(end_indices)
y_orange = total_cumulative[x_endpoints - 1]

# Fit polynomials
poly_orange = Polynomial.fit(x_endpoints, y_orange, deg=3)
print(f"Polynomial coefficients (total cumulative): {poly_orange.coef}")

# Generate smooth x values for plotting the polynomial fits
x_smooth = np.linspace(1, len(df_mono), 500)

# ============================================
# 4. Plot
# ============================================
fig, ax1 = plt.subplots(figsize=(12, 6))

# Left axis: cumulative monomial time (blue) and total cumulative (red)
ax1.fill_between(global_index + 4, total_cumulative/60, 0, color='bisque', alpha=1.0)
ax1.plot(global_index + 4, total_cumulative/60, color='darkorange', linewidth=2, label='Total cumulative time (minutes)')

ax1.fill_between(global_index + 4, cumulative_mono/60, 0,
                 color='red', alpha=1.0, label='Total cumulative time (minutes)')

ax1.set_xlabel('Global monomial index')
ax1.set_ylabel('Cumulative time (minutes)', color='black')
ax1.tick_params(axis='y', labelcolor='black')

# Polynomial fits (shifted +4 for red)
ax1.plot(x_smooth + 4, poly_orange(x_smooth)/60,
         color='darkorange', linestyle='--', linewidth=2,
         label=f'Polynomial fit (deg={poly_orange.degree()}) – total')

ax1.set_xlabel('Global monomial index')
ax1.set_ylabel('Cumulative time (minutes)', color='black')
ax1.tick_params(axis='y', labelcolor='black')

"""# Secondary axis: per‑monomial bars (orange)
ax2 = ax1.twinx()
ax2.bar(global_index + 4, df_mono['monomial_total_time_s'],
        width=0.8, color='orange', alpha=0.6, label='Time per monomial linear solve')
ax2.set_ylabel('Time per monomial (seconds)', color='orange')
ax2.tick_params(axis='y', labelcolor='orange')"""

# Vertical dashed lines at order boundaries (start of each order)
for start in start_indices:
    ax1.axvline(x=start + 4 - 0.5, color='gray', linestyle='--', alpha=0.5, linewidth=0.8)

# Combine legends
#lines1, labels1 = ax1.get_legend_handles_labels()
#lines2, labels2 = ax2.get_legend_handles_labels()
#ax1.legend(lines1 + lines2, labels1 + labels2, loc='upper left')

plt.title('Total cumulative cost = monomial steps + order steps (added at order start)')
plt.tight_layout()
plt.xlim(1, len(df_mono) + 4)
plt.ylim(0, total_cumulative.max()/60 * 1.05)
plt.grid(False)
plt.show()
# %%
