"""Plotly surface used by demo_simu.R, written to a JPEG by kaleido."""
import csv
import sys

import plotly.graph_objects as go

matrix_path, out_path, x_title, y_title, z_title = sys.argv[1:6]
with open(matrix_path, newline="") as handle:
    beta = [[float(value) for value in row] for row in csv.reader(handle)]

n_t = len(beta)
n_s = len(beta[0])
figure = go.Figure(data=[go.Surface(
    z=beta,
    x=list(range(1, n_s + 1)),
    y=list(range(1, n_t + 1)),
    showscale=False,
    colorscale=[[0, "#4575b4"], [0.5, "#fafafa"], [1, "#d73027"]],
)])
figure.update_layout(
    font=dict(size=13),
    scene=dict(
        xaxis_title=x_title,
        yaxis_title=y_title,
        zaxis_title=z_title,
        camera=dict(eye=dict(x=-1.25, y=-1.25, z=1.5)),
    ),
    margin=dict(l=0, r=0, t=0, b=0),
)
figure.write_image(out_path, width=960, height=720)
