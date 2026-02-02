# Chromosome 1 Mutation Hotspot Dashboard

This project analyses population-level variants on human chromosome 1 to:
- Calculate mutation diversity scores
- Identify mutation hotspots using sliding windows
- Extract rare variant–enriched regions
- Visualize results using an interactive dashboard

## Biological Motivation
Mutation hotspots arise due to:
- DNA replication errors
- Local sequence context
- Repair inefficiencies

These regions are important in disease genomics.

## How to Run
```bash
pip install -r requirements.txt
streamlit run app.py
