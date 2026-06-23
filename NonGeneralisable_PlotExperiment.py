# ============================================================
# Plot all experimental data from Fig. 4H sheet
# Reacted reporter signal shown as percentage
# ============================================================

from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
import xlrd

# ------------------------------------------------------------
# Load Excel sheet
# ------------------------------------------------------------
BASE_DIR = Path(__file__).resolve().parent
DATA_FILE = BASE_DIR / "data" / "in_vitro_ctRSD_modeling_data.xls"
SHEET_NAME = "Fig. 4H (SciAdv) fuel"

book = xlrd.open_workbook(DATA_FILE)
sheet = book.sheet_by_name(SHEET_NAME)

# ------------------------------------------------------------
# Data starts after the header rows
# ------------------------------------------------------------
START_ROW = 3

t_minutes = np.array(sheet.col_values(0, start_rowx=START_ROW), dtype=float)

# ------------------------------------------------------------
# Column labels
# ------------------------------------------------------------
conditions = [
    ("Input 0 nM / Fuel 0 nM", 1),
    ("Input 1.25 nM / Fuel 0 nM", 2),
    ("Input 2.5 nM / Fuel 0 nM", 3),
    ("Input 0 nM / Fuel 25 nM", 4),
    ("Input 1.25 nM / Fuel 25 nM", 5),
    ("Input 2.5 nM / Fuel 25 nM", 6),
]

# ------------------------------------------------------------
# Plot all six time courses as separate panels
# ------------------------------------------------------------
fig, axes = plt.subplots(
    2, 3,
    figsize=(12, 7),
    sharex=True,
    sharey=True
)

axes = axes.flatten()

for ax, (label, col) in zip(axes, conditions):

    y_norm = np.array(sheet.col_values(col, start_rowx=START_ROW), dtype=float)
    y_percent = 100 * y_norm

    ax.plot(t_minutes, y_percent, linewidth=2)
    ax.set_title(label)
    ax.grid(True)

# ------------------------------------------------------------
# Labels and formatting
# ------------------------------------------------------------
fig.supxlabel("Time (min)")
fig.supylabel("Reacted reporter signal (%)")
fig.suptitle("Experimental fluorescence time courses", fontsize=16)

plt.tight_layout()
plt.show()