import pandas as pd

df = pd.read_csv('hg19_chr1_variants.csv')
df.head()

import numpy as np

total = df.shape[0]
length = 248956422 #GRCh38 chr1 length
d_score = total/length
d_score

df = df.copy()
df = df.explode('alternate_bases')

transition = {('A','G'), ('G','A'), ('C','T'), ('T','C')}
def is_transition(row):
  return(row['reference_bases'], row['alternate_bases']) in transition
df['is_transition'] = df.apply(is_transition, axis = 1)

transitions = df['is_transition'].sum()
transversions = len(df) - transitions

mutation_ratio = transitions / transversions
mutation_ratio

import matplotlib.pyplot as plt

positions = df["start_position"].sort_values().values

window_size = 100
length = 248956422

window_starts = np.arange(0, length, window_size)
d_scores = []

for start in window_starts:
    end = start + window_size
    count = np.sum((positions >= start) & (positions < end))
    d_scores.append(count / window_size)

df_windows = pd.DataFrame({
    "window_start": window_starts,
    "diversity": d_scores
})

mean = df_windows["diversity"].mean()
std = df_windows["diversity"].std()

hotspot_threshold = mean + 2 * std

df_windows["hotspot"] = df_windows["diversity"] > hotspot_threshold

plt.figure(figsize=(15, 5))

plt.plot(
    df_windows["window_start"],
    df_windows["diversity"],
    linewidth=1,
    label="Diversity"
)

plt.scatter(
    df_windows.loc[df_windows["hotspot"], "window_start"],
    df_windows.loc[df_windows["hotspot"], "diversity"],
    label="Hotspots"
)

plt.xlabel("Chromosome 1 position (bp)")
plt.ylabel("Diversity score (variants / bp)")
plt.title("100-bp Diversity on Chromosome 1")
plt.legend()

plt.show()

rare_threshold = 0.01

df_rare = df[df["MAF"] < rare_threshold].copy()
df_rare.shape

window_size = 100

df_rare["window"] = (df_rare["start_position"] // window_size) * window_size

rare_window_counts = (
    df_rare
    .groupby("window")
    .size()
    .reset_index(name="rare_variant_count")
)

rare_window_counts["mutation_rate"] = (
    rare_window_counts["rare_variant_count"] / window_size
)

above_avg_windows = rare_window_counts[
    rare_window_counts["mutation_rate"] > d_score
]

above_avg_windows.head()

import streamlit as st
import pandas as pd
import matplotlib.pyplot as plt

st.title("Chromosome 1 Mutation Hotspot Dashboard")

st.markdown("Public population genomics – 100 bp resolution")

threshold = st.slider(
    "Hotspot threshold (mean + N×SD)",
    1.0, 4.0, 2.0, 0.5
)

cutoff = df_windows["diversity"].mean() + threshold * df_windows["diversity"].std()
df_windows["hotspot"] = df_windows["diversity"] > cutoff

fig, ax = plt.subplots(figsize=(12,4))
ax.plot(df_windows["window_start"], df_windows["diversity"])
ax.scatter(
    df_windows[df_windows["hotspot"]]["window_start"],
    df_windows[df_windows["hotspot"]]["diversity"]
)
ax.set_xlabel("Chr1 position (bp)")
ax.set_ylabel("Diversity")
st.pyplot(fig)

st.subheader("Rare Variant Hotspot Windows")
st.dataframe(above_avg_windows.sort_values("mutation_rate", ascending=False))
