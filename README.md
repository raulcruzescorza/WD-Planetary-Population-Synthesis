# Population Synthesis of Planetary Systems around White Dwarfs

> This repository contains the Python code for the MPhil thesis "Surviving Worlds: A Stochastic Model of Planet Demographics around Future White Dwarfs". The project consists of a self-consistent population synthesis simulation designed to model the formation and evolution of planetary systems in the solar neighborhood, with the goals of predicting the final demographics of planets (with $\leq4$ R🜨) that survive around white dwarfs and establishing a fraction of these stars with at least one planet.

<p align="center">
  <img src="path/to/your/2d_occurrence_map.png" width="700" alt="2D Occurrence Map of Surviving Planets">
</p>

## Overview of the Model

The simulation follows a multi-step process to generate and evolve a realistic population of stars and planets:

* **Star Formation History (SFH):** A stellar population is generated over a 14-Gyr history. The model uses the cosmic star formation rate from **(Madau & Dickinson, 2014)**, normalised to match the observed local stellar census of ~5,230 stars within 25 pc from **(Golovin et al., 2023)**.

* **Stellar Population:** Stars are assigned masses from a composite Chabrier/Salpeter Initial Mass Function (IMF). Although the IMF generates a wide mass spectrum (0.08-8M⊙), the simulation focuses on 0.8–2.0 M⊙ progenitors. Lifetimes are calculated using the formulae from **(Hurley et al., 2000)**.

* **Planetary Population:** Each progenitor star is populated with planets based on the occurrence rates from **(Hsu et al., 2019)**. The number of planets per star is determined stochastically from a Poisson distribution whose mean is the total occurrence rate.

* **Planetary System Evolution:** As stars reach the end of their main-sequence lifetime, the orbits of their planets are evolved numerically with `scipy.integrate.solve_ivp`. The model accounts for the two competing effects:
    * **Stellar Mass Loss:** Causes orbits to expand.
    * **Tidal Forces:** Cause orbits to shrink and can lead to engulfment.
    The implementation is based on the model in **(Sanderson et al., 2022)**, and uses the stellar evolution tracks from **(Vassiliadis & Wood, 1993)**.

## Key Results

The main scientific outputs of this simulation are:

* A validation of the model by successfully reproducing the observed number of white dwarfs in the 25 an 40 pc solar neighborhood.
* A full census of the surviving planet population, including their final orbital distances and radii.
* A 2D occurrence map showing the probability of finding a surviving planet as a function of its mass and orbit.
* An optimistic upper limit on the fraction of white dwarfs that could host planets in their habitable zones.

## Getting Started

Follow these steps to set up and run the code on your local machine.

### Prerequisites

* Python 3.x
* A `requirements.txt` file is included in the repository to install all necessary dependencies.


### Installation

1.  **Clone this repository:**
    ```bash
    git clone https://github.com/raulcruzescorza/WD-Planetary-Population-Synthesis.git
    cd WD-Planetary-Population-Synthesis
    ```

2.  **Navigate to the project directory**
    ```bash
    cd WD-Planetary-Population-Synthesis
    ```

3.  **Create a virtual environment**

    This command creates a directory named `venv` which will contain a private Python installation and the necessary libraries for this project.
    ```bash
    python -m venv venv
    ```

4.  **Activate the virtual environment**

    You need to activate the environment in every new terminal session you open to work on this project.

    * **On macOS / Linux:**
        ```bash
        source venv/bin/activate
        ```
    * **On Windows:**
        ```bash
        .\venv\Scripts\activate
        ```
    > You'll know it worked because your terminal prompt will change to show `(venv)` at the beginning.

5.  **Install the required dependencies**

    This command reads the `requirements.txt` file and automatically installs all the specific library versions used for this analysis.
    ```bash
    pip install -r requirements.txt
    ```

6.  **Run the project**

    You are all set! Now you can launch Jupyter and run the `FinalNotebook.ipynb` file.
    ```bash
    jupyter notebook
    ```
    Or if you use Jupyter Lab:
    ```bash
    jupyter lab
    ```
7.  **End of the virtual environment**

    To return to your kernel settings use:
    ```bash
    deactivate
    ```
    Only use this line when you are done running everything.

### Running the Simulation

1.  Ensure all input data files (e.g., `agb1p0.dat`, `ajab31abt2_ascii.txt`) are in the correct directory.
2.  Launch Jupyter Notebook or Jupyter Lab:
    ```bash
    jupyter notebook
    ```
3.  Open and run the cells in `FinalNotebook.ipynb` sequentially. The notebook is divided into sections for data processing, simulation, and analysis.

## Citation

If you use this code or its results in your research, please cite the original thesis:

* R. A. Cruz Escorza (2025), *Surviving Worlds: A Stochastic Model of Planet Demographics around Future White Dwarfs*, MPhil Thesis, University of Cambridge.

## License

This project is distributed under the MIT License. See the `LICENSE` file for more details.

## Author

* **R. A. Cruz Escorza** - [GitHub](https://github.com/your-username) - [Email](rac228@cam.ac.uk)
