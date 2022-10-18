# Bottom baryons decay widths and mass spectra

Code to compute bottom baryon spectra and decay widths. A fit is performed to obtain the model parameters. Errors are propagate via bootstrap Monte Carlo Gaussian sampling

## Framework installation

To install the framework you need anaconda and git on a linux machine. In the terminal(powershell):
1. Clone the repository:
  ```
  git clone git@github.com:andrex-naranjas/bottom-baryonsFW.git
  ```
2. Access the code:
  ```
  cd bottom-baryonsFW
  ```
3. Install the conda enviroment:
  ```
  conda env create -f config.yml
  conda activate bottom
  conda develop .
  ```
4. Compile the decay widths C++ code (here we use C++17):
  ```
  cd ./decays/DecayWidths/
  make obj
  make
  ```
5. Minimal run:
  ```
  python3 ./scripts/bootstrap_three_quark.py omegas
  python3 ./scripts/bootstrap_diquark.py omegas
  python3 ./scripts/print_results.py omegas
  ```
X. Check that your plots ./plots directory