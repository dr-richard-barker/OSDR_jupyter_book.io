---
description: >-
  The OSDR RadLab application is a valuable tool for researchers and scientists
  who need access to radiation measurements from space, including astrobotanists
  interpreting how plants respond to the spaceflight radiation environment.
---

# RadLab Overview

**Please follow this link to access the** [**RadLab OSDR Visualization application**](https://visualization.osdr.nasa.gov/radlab/gui/data-overview/)**.**

{% embed url="https://visualization.osdr.nasa.gov/radlab/gui/data-overview/" %}

The application provides a central repository for radiation data from a wide variety of sensors that have flown on various spacecraft. This makes it easy for users to find and access the data they need, regardless of the specific sensor or spacecraft that collected it.

The OSDR RadLab application is also a powerful tool for data analysis. It provides a variety of features that allow users to visualize, analyze, and compare radiation data. This can be helpful for identifying trends and patterns in the data, as well as for understanding the effects of radiation on different materials and systems. Ionizing radiation is one of the defining environmental factors of spaceflight, and it affects all living organisms — including plants. For astrobotany researchers, the dose and timing of radiation a payload experienced is valuable context for interpreting how plants such as *Arabidopsis thaliana* responded during a mission.

One of the most important features of the OSDR RadLab application is its ability to connect primary radiation data with metadata such as vehicle and time. This information is essential for understanding the context in which the radiation data was collected. It can also be used to filter and sort the data, making it easier to find the specific information that users are looking for.

The OSDR RadLab application is a valuable resource for researchers and scientists who need access to radiation measurements from space. It is a powerful tool for data analysis and visualization, and it can help users to understand the effects of radiation on different materials and systems. Plant space-biology investigations benefit from this too: spaceflight experiments such as the Characterizing Arabidopsis Root Attractions (CARA) study (OSD-120) examine how *Arabidopsis* grows in orbit, and pairing those biological results with RadLab's radiation records helps researchers account for the radiation environment the plants were exposed to.

#### **Data overview page** <a href="#txgh0vy0ejzk" id="txgh0vy0ejzk"></a>

**When you first land on the RadLab you start in the “data overview” page.**

This page summarizes the different sensors, the time they were collecting data and the space vehical they were on. The graph is interactive so the users can select subsets within the data to drill down and the graph will rearrange to optimize the data visualization. On the bar on the left the user can then navigate down to a range of different visualization options.

#### RadLab Time series line plots <a href="#j9p2w6eomst7" id="j9p2w6eomst7"></a>

**When you first land on the Time series plot page the graphs initially appear empty**

This page allows users to view the radiation data but requires users needs to define the periods and sensors of interest.

* In this example, a Red arrow highlights how users can quickly load a preselected period to help get familiarized with the interface.
* Select either “Total dose rate” or “Total flux” values.
* The user can alternate between a “Linear” and “Log” scale.


---

## Embedded RadLab notebook with generated plots

To make the RadLab overview more useful, this book includes a companion executed notebook that pulls live RadLab telemetry and renders interactive Plotly graphs (time series, mission comparisons, cumulative dose, and Earth vs space baselines).

- Notebook (source + executed outputs): `chapters/OSDR_radlab_astrobotany.ipynb`
- Built book page: `OSDR_radlab_astrobotany.html` (available in the built site alongside this chapter)

View the interactive notebook in the book:

[Open the RadLab radiation notebook — "Radiation environment of OSDR astrobotany missions"](OSDR_radlab_astrobotany.html)

You can also embed the notebook page here (if your site allows iframes):

<iframe src="/OSDR_radlab_astrobotany.html" width="100%" height="800" style="border:1px solid #ddd"></iframe>

Notes for maintainers and contributors

1. To regenerate the notebook locally (it pulls live telemetry from the RadLab API):

```bash
python chapters/_make_radlab_notebook.py
jupyter nbconvert --to notebook --execute chapters/OSDR_radlab_astrobotany.ipynb
```

2. The notebook uses Plotly for interactive figures. When executed in the jupyter-book build environment you may need to ensure the renderer is appropriate for your build environment (the notebook currently sets `pio.renderers.default = 'notebook_connected'`).

3. If you prefer static PNG/SVG figures for CI or static builds, execute the notebook and save static images from the Plotly objects (the notebook includes commented examples for exporting static figures).

If you'd like, I can (a) execute the notebook in the repository and commit the updated notebook with outputs, or (b) modify this page to include a few inline static example images generated from RadLab telemetry. Which would you prefer?