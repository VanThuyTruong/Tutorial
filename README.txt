README

About the Project
This project is simulates tumour cell survival during targeted 
anti-cancer therapy with the MEK inhibitor cobimetinib using 
an agent-based model and an ordinary differential equation model.


Built With
Python 3.7


Prerequisites
Installed python

Required python packages:
numpy
os 
random
matplotlib
scipy
time
math
PIL
opencv
PPP

Getting Started
To start the simulation run simulation.py.

For changing parameters and models use config.py.

abm.py runs the agent-based model.
abm_plots.py creates the agent-based model plots.

ode.py runs the ordenary differential equation model.
ode_plots.py creates the ordenary differential equation model.

cell.py sets up a tumour cell as one agent.
cells.py sets up a tumour cell population consisting of tumour cells.
with different pERK values as attributes, in addition it codes the 
behaviour rules of the tumour cells.

All the python files and the PPP package folder need to be stored in the 
same folder for running the code.

Usage
Simulation of different dosage regimen (multiple/single dose, different 
drug amounts given ect.), different pharmacokinetic and pharmacodynamic 
parameters, different pERK distributions (uniform, bimodal, normal ect.), 
comparison between model behaviour of the agent-based model and ordinary 
differential equation model


License
Creative Commons Attribution (CC-BY) license 

Contact
Van Thuy Truong
PhD Researcher in Clinical Pharmacology Biologics and Bioanalysis 
vn.thuy.truong@gmail.com


Acknowledgements
This project is part of the QuanTII (Quantitative T cell Immunology 
and Immunology) project funded by the Horizon 2020 programme of the 
European Union.

Thanks to Grant Lythe and Joseph Elmes for helping to optimize the code.