"""Plot the compact double-mutant plan summary."""
from __future__ import annotations
import csv
from collections import Counter
from pathlib import Path
from typing import Mapping,Sequence
import matplotlib.pyplot as plt
def plot_plan(rows:Sequence[Mapping[str,object]],data_path:Path,png_path:Path,svg_path:Path)->None:
 counts=Counter(str(row["combination_track"]) for row in rows);data=[{"combination_track":key,"candidate_count":counts[key]} for key in ("affinity_x_affinity","affinity_x_property","property_x_property")]
 with data_path.open("w",encoding="utf-8-sig",newline="") as h:w=csv.DictWriter(h,fieldnames=list(data[0]),lineterminator="\n");w.writeheader();w.writerows(data)
 labels=["Affinity × affinity","Affinity × property","Property × property"];values=[row["candidate_count"] for row in data]
 fig,ax=plt.subplots(figsize=(6.4,4.1));bars=ax.bar(labels,values,color=["#4477AA","#228833","#CCBB44"]);ax.set_ylabel("Unfiltered double-mutant count");ax.set_ylim(0,max(values)*1.18);ax.spines[["top","right"]].set_visible(False)
 for bar,value in zip(bars,values,strict=True):ax.text(bar.get_x()+bar.get_width()/2,value+1,str(value),ha="center",va="bottom")
 ax.text(0.01,0.98,"86 valid pairs retained; 5 same-position pairs excluded",transform=ax.transAxes,va="top",fontsize=9);fig.tight_layout();fig.savefig(png_path,dpi=600);fig.savefig(svg_path);plt.close(fig)
