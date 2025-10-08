# VeBaMoS: Vector-Based Molecular Similarity

**Authors:** Camilo Prada Latorre, Ricardo Vivas-Reyes, and Jhon E. Zapata-Rivera  
**Affiliations:** Universidad de los Andes · Universidad de Cartagena · Universidad del Valle  

---

## 🧠 Overview

**VeBaMoS (Vector-Based Molecular Similarity)** is a Python-based framework that quantifies molecular similarity by representing molecules as descriptor-based vectors and comparing them through a unified similarity index.  

Unlike many existing methods that rely only on geometric or electronic properties, **VeBaMoS integrates both**, yielding a chemically meaningful, explainable, and flexible definition of molecular similarity.  

> 🧪 Developed as part of academic research at the Universidad de los Andes and Universidad del Valle.

---

## ✨ Features

- **Hybrid descriptor space:** combines geometric and electronic molecular information.  
- **Vector-based formalism:** interpretable similarity in terms of magnitude and direction.  
- **Normalized index:** values between 0 (dissimilar) and 1 (identical).  
- **Quantum integration:** reads data directly from ORCA output files.  
- **Command-line interface (CLI):** flexible and scriptable.  
- **Structured output:** results stored as `.smo` files.  
- **Customizable descriptors:** users select which properties to include.  

---

## ⚙️ Mathematical Foundations

Each molecule is represented as a vector in a descriptor space:

\[
\vec{X} = \sum_{i=1}^{n} x_i \vec{a_i}
\]

To address the difference in magnitude between each descriptor, they are normalized by converting it into a fraction of the total sum of its values. After normalization:

\[
x'_i = \frac{x_i}{\sum_j x_j}
\]

The similarity between two molecules \( X \) and \( Y \) is quantified by:

\[
\begin{aligned}
\alpha &= ||\vec{X}|| - ||\vec{Y}|| \\
\beta &= \arccos\left(\frac{\vec{X}\cdot\vec{Y}}{||\vec{X}|| \, ||\vec{Y}||}\right)
\end{aligned}
\]

We then calculate a final similarity index to improve the interpretability:

\[
\sigma(z) = \frac{1}{e^{z / \log_{10}(N)}}, \quad z = |\alpha| + |\beta|
\]

where \( N \) is the number of molecules in the dataset.

---

## 🧩 Installation

### Requirements

- **OS:** Linux (recommended)
- **Python:** ≥ 3.8  
- **Quantum package:** [ORCA 4.2.1+](https://orcaforum.kofo.mpg.de)  
- **Dependencies:**
  ```bash
  numpy
  scipy
  pandas
  matplotlib
  os
  sys
  ```

### Installation Steps

1. Clone this repository:
   ```bash
   git clone https://github.com/CPradaL/VeBaMoS.git
   cd VeBaMoS
   ```

2. Install dependencies:
   ```bash
   pip install numpy
   pip install scipy
   pip install matplotlib
   pip install pandas
   ```

3. Add a shell alias (optional, for convenience):
   ```bash
   echo "alias vebamos='python ~/VeBaMoS/VeBaMoS.py'" >> ~/.bashrc
   source ~/.bashrc
   ```

4. Run VeBaMoS:
   ```bash
   vebamos -i <input_file>
   ```

---

## 📥 Input Format

VeBaMoS uses a **structured text input file** to specify settings, molecules, and descriptors.

Example (`input.inp`):

```bash
!Out=(Results) Vectores=(Norm) ExpInf=(activities.txt,testosterone)
!Inp=(/home/user/molecules)

%archivos
mol1
mol2
mol3
end

%Descriptores
dipole
volume
HOMO-LUMO
end
```

### Keywords
| Keyword | Description |
|----------|--------------|
| `Out` | Defines output filename (saves as `.smo`) |
| `Vectores` | Type of vector printed (`Full` or `Norm`) |
| `Inp` | Path to directory containing molecular files |

### Blocks
| Block | Purpose |
|--------|----------|
| `%archivos ... end` | List of molecule files to compare (e.g. `molX.out` ) |
| `%Descriptores ... end` | Optional: specify which descriptors to include |

Each molecule must have a corresponding **ORCA output file** (`molX.out`) containing geometry and electronic properties.

---

## 🧮 Available Descriptors


| Input File Command | Descriptor |
|--------|----------|
| Dipolarx | Dipole moment (x axis) |
| Dipolary | Dipole moment (y axis) |
| Dipolarz | Dipole moment (z axis) |
| DipolarT | Total Dipole moment |
| E_Disper | Dispersion Energy |
| Distanc1 | Longest interatomic distance |
| VolMolec | Molecular Volume |
| A_Superf | Superficial Area |
| Qdrupolx | Quadrupole moment (x axis) |
| Qdrupoly | Quadrupole moment (y axis) |
| Qdrupolz | Quadrupole moment (z axis) |
| QdrupolI | Isotropic Quadrupole moment |
| Polarizx | Polarizability (x axis) |
| Polarizy | Polarizability (y axis) |
| Polarizz | Polarizability (z axis) |
| PolarizI | Isotropic Polarizability |
| EnerHOMO | HOMO Energy |
| EnerLUMO | LUMO Energy |
| Distanc2 | Longest distance from any atom to the line defined by D1 |
| Gap_HOLU | HOMO-LUMO Gap |
| Dureza__ | Chemical Hardness |
| Suavidad | Chemical Softness |
| PQuimico | Chemical Potential |
| Elfilici | Global Electrophilicity |


-  
- Quadrupole moment  
- Polarizability  
- Molecular volume  
- Surface area  
- Maximum molecular length  
- Perpendicular molecular length  
- Ionization energy  
- Electron affinity  
- HOMO–LUMO gap  
- Global hardness / softness  
- Chemical potential  
- Electrophilicity  

---

## ⚖️ License

VeBaMoS is distributed for **academic and research use only**.  
Please contact the authors before redistribution or modification.

---

## 📫 Contact

**Main Developer:**  
Camilo Prada Latorre — [c.prada@uniandes.edu.co](mailto:c.prada@uniandes.edu.co)  

**Advisor:**  
Jhon E. Zapata-Rivera — [jhon.zapata.rivera@correounivalle.edu.co](mailto:jhon.zapata.rivera@correounivalle.edu.co)
