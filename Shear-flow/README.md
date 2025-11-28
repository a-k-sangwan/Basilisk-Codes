# 🌊 Two-Phase Shear Flow Simulation using Basilisk

This repository contains scripts and simulation code to study two-phase shear flow instability using the **Basilisk** CFD framework. The simulation initializes a perturbed air–water interface and tracks its evolution under shear and gravity effects.

## 📂 Repository Structure

```
single.c          # Main Basilisk simulation
output_vtu.h      # VTU output utilities for ParaView visualization
genData.py        # Initial velocity & interface generator
run.sh                 # Compile & run wrapper script
.gitignore
LICENSE
README.md
```

> Output folders (`dumpfile/`, `vtufile/`, `interface/`) are automatically generated during simulation.

---

## 🔧 Requirements

### Basilisk
You must install Basilisk and its dependencies.  
Instructions: https://basilisk.fr/src/INSTALL

Core Basilisk modules used:
- `navier-stokes/centered.h`
- `two-phase.h`
- `reduced.h`
- `tension.h`
- `distance.h`
- `output_vtu_foreach.h`

### Python Dependencies
Run:
```bash
pip3 install numpy sympy scipy matplotlib
```

---

## ▶️ How to Run

### Option 1 — Use provided script
```bash
chmod +x run.sh
./run.sh
```

### Option 2 — Manual compile & run
```bash
qcc -O2 -Wall single.c -lm -o run
./run
```

Simulation automatically:
✔ Creates output folders  
✔ Generates initial fields  
✔ Exports VTU snapshots  
✔ Logs runtime info

---

## 🧪 Physics Overview

Simulation setup:
- **Two-phase incompressible Navier–Stokes**
- **Air over water** with **density ratio = 0.5**
- **Surface tension = 72**
- **Sheared top boundary** (`u.t[top] = U0`)
- **Gravity downwards (`G.y = -980`)**
- **Periodic domain horizontally**
- **No viscosity** for both fluids

Grid:
- Initial resolution: 512
- Adaptive refinement up to `MAXLEVEL = 10`

---

## 📤 Output Files

| File Type | Location | Purpose |
|----------|----------|---------|
| VTU snapshots | `vtufile/` | Velocity + pressure + scalar fields for ParaView |
| Interface geometry | `interface/` | Visualize VOF interface |
| Restart dumps | `dumpfile/` | Used to continue simulations |
| Timelog | `track.out` | Diagnostics |

To view in **ParaView**:
```
Open > vtufile/snap-*.vtu
```

---

## 🔁 Restarting a Simulation

Say you want to restart the simulation from 20th dumpfile i.e., File dumpfile/dump-20 
```bash
./run 20
```

---

## 📈 Python Field Generation

`python/genData.py` computes:
- Initial shear velocity distribution
- Perturbed interfacial profile η(x)

Usage:
```bash
python3 genData.py -vel   # Generate initial velocity
python3 genData.py -eta   # Generate interface shape
```

---

## ⚠️ Notes

- This configuration is set up for **2D**
- VTU & dump intervals controlled by: `tdump` variable
- Simulation final time controlled by: `t_final`
- Memory and runtime scale with grid & refinement

---

## 📜 License

MIT License — free to use, modify, and distribute.

---

## 👨‍💻 Author

Anil Kumar, Chemical Engineering, Indian Institute of Technology Bombay, India 400 076

---

## How to cite / acknowledge

If you use this repository in research or publications, please acknowledge:

- The Basilisk C framework (https://basilisk.fr)

- This GitHub repository

### Happy Simulating with Basilisk! 🌊🚀
