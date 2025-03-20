## Workflow for Cross-Section Extraction

### Preselection from HIPO File

This step is performed by `eIDcut_DCfiducial.cpp` in my workflow. It reads through all event banks and stores the events that pass the preselection cut at lines `875–886`.

The output of this file is a ROOT file containing all the necessary variables.

To run this code:

```
source go
./eIDcut_DCfiducial
```

This preselection code is typically what I consider the `main.cpp`. Ideally, we apply all necessary cuts loosely at this stage to minimize the amount of data that needs to be stored while avoiding the loss of edge values that may be useful later. For example, if the optimal vertex cut is \( z > -10 \), we might apply \( z > -14 \) here. The subsequent analysis will be based on the data stored from this step.

### Analysis Code

An example of an analysis script is `Draw_Electron_Kinematics.cpp`. This type of code usually includes predefined functions for cuts and corrections before performing the analysis and generating histograms.

To run this code:

```
clas12root Draw_Electron_Kinematics.cpp
```

However, I believe it’s better to define all cuts and corrections in separate function files that can be called by every analysis script. This makes them more reusable and easier to manage. Additionally, when writing functions for cuts and corrections, I recommend setting different cut/correction parameters as function variables instead of hardcoding them. This way, we don’t need to modify the function files directly when changing parameters. For example, in the following function:

$$
\Delta P = \text{para}[1] \cdot \theta^2 + \text{para}[2] \cdot \theta + \text{para}[3]
$$

We can define it as:

```
proton_momentum_correction(event, para[1], para[2], para[3])
```

### Binning Code

To create bins for different kinematic ranges and extract the yield, I use `auto_bin_analysis_W_corrected.cpp`.

This script processes the electron yield separately for different sectors and $\theta$ bins. We should always incorporate this binning step as part of our workflow. I believe this step is crucial, though I forgot to mention it during our meeting.

To run this code:

```
clas12root auto_bin_analysis_W_corrected.cpp
```

### Normalization Code

After obtaining the yield for different bins, the next step is normalization.

For normalization, we use `normalization.cpp`. Although normalization could be incorporated into the analysis or binning code, I prefer to run it separately. This is because normalization reads histograms from previous steps rather than raw data from the preselection step. Since reading the full dataset can be time-consuming, running normalization separately allows us to modify output styles or formats without having to reprocess everything from scratch.

To run this code:

```
clas12root normalization.cpp
```

I believe these four steps are the most crucial parts of the workflow. Once we complete these steps, the rest should be straightforward.
