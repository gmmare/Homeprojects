import plotly.graph_objects as go

fig = go.Figure(data=[go.Sankey(
    node = dict(
      pad = 15,
      thickness = 20,
      line = dict(color = "black", width = 0.5),
        label=["", "", "", "", "", "",
               "", ""],  # 8 items
      color = "blue"
    ),
    link = dict(
      source = [0,0,1,1,2,2,3,3,4,4], # indices correspond to labels, eg A1, A2, A2, B1, ...
      target = [5,1,2,6,3,6,4,7,5,8],
      value = [1,99,60.46,(99 - 60.46),55.24,(60.46 - 55.24),39.01,(55.24 - 39.01)]
  ))])

fig.update_layout(title_text="Basic Sankey Diagram", font_size=10)
fig.show()

# label = ["Chemical Power", "Heat", "Gas Power", "Propulsive Jet Power", "Thrust Power",
#           "Incomplete combustion",
#           "Heat","Kinetic Energy"], #8 items