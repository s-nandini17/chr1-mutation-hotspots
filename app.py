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
