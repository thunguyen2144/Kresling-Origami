# Kresling Origami Simulator

A Python/JAX implementation for simulating Kresling origami structures using finite element analysis.

![Kresling Structure](results/kresling.jpg)

---

## Quick Start

### 1. Setup Environment

**Create conda environment:**
```bash
conda env create -f environment.yml
```

This will create an environment named `kresling` with:
- Python 3.10
- JAX (CPU-only, works on Windows/Linux/macOS)
- NumPy 1.x (for compatibility)
- Matplotlib, Pillow, SciPy

**Activate the environment:**
```bash
conda activate kresling
```

**Note:** You must activate the environment every time you want to run simulations.

### 2. Run Simulation

```bash
cd Desktop\Kresling-Simulator-main
cd Desktop\test_data
conda activate kresling
python main.py
```

- **Solver:** RK4 (explicit, 4th-order accuracy) - **Default solver**
- **Features:** Time-varying sinusoidal loads, Kelvin-Voigt damping
- **Runtime:** 5-10 minutes
- **Output:** All results saved to `results/` folder

**Note:** First run is slower due to JAX JIT compilation

### 3. View Results

Check the `results/` folder for generated files:
- `kresling_final_shape.png` - Final deformed structure
- `kresling_velocity.png` - Velocity time history  
- `kresling_displacement.png` - Displacement history
- `kresling_force_history.png` - Force time history
- `Animation.gif` - Animated deformation with force vectors

---

## Adjusting Parameters

Edit `main.py` to customize your simulation:

### Geometry Parameters (line ~50-54)
```python
R = 40e-3           # Radius: 40 mm → try 30-100 mm
H = 50e-3           # Height: 50 mm → try 30-150 mm
theta = 30 / 180 * np.pi  # Twist: 30° → try 15-45°
N = 6               # Sides: 6 (hexagon) → try 3-12
m = 1               # Stories: 1 → try 1-5
```

### Material Parameters (line ~57-61)
```python
spr_stiff = 1e-7    # Spring stiffness (N·m)
bar_E = 1e7         # Young's modulus (Pa)
density = 700       # Density (kg/m³)
panel_v = 0.3       # Poisson's ratio
panel_t = 2e-3      # Thickness: 2 mm
```

### Damping Parameters (line ~64-65)
```python
bar_eta = 1e3       # Bar viscosity (Pa·s)
spr_c_rot = 1e-4    # Spring rotational damping (N·m·s)
```

### Time-Varying Loading (line ~109-111)
```python
force_amplitude = 20.0      # Force amplitude (N)
omega = 200 * np.pi         # Angular frequency (rad/s)
load_frequency = 100        # Hz
```

### Simulation Settings (line ~134-135)
```python
dt = 1e-5           # Time step (seconds)
num_steps = 10000   # Number of steps
```

**See [PARAMETERS_GUIDE.md](PARAMETERS_GUIDE.md) for complete parameter reference.**

---

## Example: Running the Simulation

```bash
conda activate kresling
python main.py
```

**Features:**
- RK4 solver (explicit, 4th-order accuracy) - **Default**
- Time-varying sinusoidal loading (customizable frequency)
- Kelvin-Voigt material damping model
- Force vector visualization in animation
- All results saved to `results/` folder

**Note:** RK4 is the default solver. Newmark-beta is available in the codebase for development and comparison purposes.

**Output files:**
- `results/kresling_final_shape.png` - Final deformed structure
- `results/kresling_velocity.png` - Velocity time history
- `results/kresling_displacement.png` - Displacement history
- `results/kresling_force_history.png` - Force time history
- `results/Animation.gif` - Animation with force vectors

### Example Parameter Modifications

#### Octagonal Structure
```python
N = 8               # 8 sides instead of 6
R = 60e-3           # 60 mm radius
```

#### Multi-Story Tower
```python
m = 3               # 3 stories instead of 1
H = 70e-3           # 70 mm per story
```

#### Higher Frequency Loading
```python
force_amplitude = 50.0    # Increase force amplitude
omega = 500 * np.pi       # Higher frequency
```

---

## Project Structure

```
kresling_simulator/
├── README.md                # This file
├── PARAMETERS_GUIDE.md      # Complete parameter reference
├── environment.yml          # Conda environment specification
├── main.py                  # Main simulation script
│
├── kresling/                # Main package
│   ├── geometry/            # Geometry generation
│   ├── fem/                 # Finite element analysis
│   ├── solvers/             # Time integration
│   └── visualization/       # Plotting and animation
│
└── results/                 # Output files (generated)
```

---

## Prerequisites

- **Conda** (Miniconda or Anaconda) - Required
  - Download: https://docs.conda.io/en/latest/miniconda.html
  - Install it before proceeding

---

## Alternative: Manual Installation (without environment.yml)

If you prefer not to use the environment.yml file:

```bash
# Create environment manually
conda create -n kresling python=3.10
conda activate kresling

# Install dependencies with correct NumPy version
pip install "jax[cpu]" "numpy>=1.24.0,<2.0.0" matplotlib pillow scipy
```

**Important Notes:**
- **NumPy version:** Must be `<2.0.0` for compatibility with matplotlib
- **JAX on Windows:** Must use `jax[cpu]` (no GPU support on Windows)
- **JAX on Linux/macOS:** Can also use `jax[cpu]` for CPU-only execution

---

## Troubleshooting

### Conda not found
**Problem:** `conda: command not found`  
**Solution:** Install Conda from https://docs.conda.io/en/latest/miniconda.html and restart terminal

### Wrong environment active
**Problem:** Running from wrong environment (e.g., `bldc` instead of `kresling`)  
**Solution:** 
```bash
conda activate kresling
```
Check your prompt - it should show `(kresling)` at the beginning

### NumPy compatibility error
**Problem:** "A module that was compiled using NumPy 1.x cannot be run in NumPy 2.x"  
**Solution:** 
```bash
# Option 1: Recreate environment (recommended)
conda env remove -n kresling
conda env create -f environment.yml
conda activate kresling

# Option 2: Fix current environment
conda activate kresling
pip install "numpy>=1.24.0,<2.0.0" --force-reinstall
```

### Environment already exists
**Problem:** `CondaValueError: prefix already exists`  
**Solution:**
```bash
# Remove old environment
conda env remove -n kresling

# Create fresh environment
conda env create -f environment.yml
```

### Simulation is slow
**This is normal:** First run is slower (JAX compiling functions)  
**For testing:** Set `num_steps = 100` in `main.py`

### Import errors
**Problem:** `ModuleNotFoundError` or import failures  
**Solution:** Make sure environment is activated:
```bash
conda activate kresling
python -c "import jax, numpy, matplotlib; print('All imports OK')"
```

---

## Features

- 🚀 **JAX acceleration** - JIT compilation for performance
- 📊 **Rich visualization** - 3D plots and GIF animations
- 🔧 **Easy customization** - Clear parameter sections
- 🪟 **Cross-platform** - Windows, Linux, macOS
- 🐍 **Clean code** - Modular, well-documented

---

## Technical Details

### Finite Element Analysis
- **Bar elements** - Linear elastic with geometric nonlinearity
- **Rotational springs** - 4-node dihedral angle springs
- **Damping** - Kelvin-Voigt material damping

### Time Integration Solvers

#### RK4 (Default)
The 4th-order Runge-Kutta method is the **default solver** used in `main.py`:
- **Type:** Explicit time integration
- **Accuracy:** 4th-order (O(Δt⁴))
- **Advantages:** 
  - High accuracy for dynamic problems
  - Efficient for systems with moderate nonlinearity
  - Well-suited for time-varying loads
- **Best for:** General-purpose simulations, time-varying loading, dynamic analysis

#### Newmark-beta (Optional)
Available in the codebase for development and comparison purposes:
- **Type:** Implicit time integration
- **Accuracy:** 2nd-order (O(Δt²))
- **Advantages:**
  - Unconditionally stable (β=1/4, γ=1/2)
  - Better for stiff systems
- **Use case:** Development, testing, comparison studies

**To switch solvers:** Import `solve_newmark_beta` instead of `solve_runge_kutta_4` in `main.py`

### Mathematical Model
Dynamic equilibrium: **M**ü + **C**u̇ + **K**u = **F**(t)

Where:
- **M**: Mass matrix (lumped)
- **C**: Damping matrix (Kelvin-Voigt)
- **K**: Stiffness matrix (material + geometric)
- **F**: External force vector

---

## Performance

**Typical runtimes on modern CPU:**
- **Quick test** (num_steps=1000): ~2 minutes
- **Standard** (num_steps=10000): 5-10 minutes
- **Extended** (num_steps=20000): 15-20 minutes

**Features:**
- Includes time-varying force visualization
- Multiple plots (displacement, velocity, force history)
- Animated GIF with force vectors
- All results saved to `results/` folder

**Note:** First run includes JAX JIT compilation overhead (~30 seconds).

---

## Documentation

- **README.md** (this file) - Quick start and overview
- **PARAMETERS_GUIDE.md** - Complete parameter reference with examples
- **CHANGES_SUMMARY.md** - Technical implementation details
- **DAMPING_PARAMETERS.md** - Damping model documentation

---

## Platform Support

- ✅ **Windows** - CPU-only JAX (fully functional)
- ✅ **Linux** - CPU JAX (fully functional)
- ✅ **macOS** - CPU JAX (optimized for Apple Silicon)

---

## License

[Specify your license]

---

## Citation

If you use this code in research, please cite:
```
[Add citation information]
```

---

## Contact

[Your contact information]

---

**Ready to start?**

```bash
# 1. Create environment
conda env create -f environment.yml

# 2. Activate it
conda activate kresling

# 3. Run simulation
python main.py

# 4. Check results/ folder for outputs!
```

