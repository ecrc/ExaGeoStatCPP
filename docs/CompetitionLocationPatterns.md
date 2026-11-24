# Competition-inspired location patterns

Synthetic generation accepts `--location-pattern=jittered-grid|regular|chessboard|left-bottom|clusters`. The default retains the existing jittered-grid generator. Non-default patterns support spatial 2D only; covariance kernel and theta remain separate settings. Both StarPU generation and PaRSEC's actual STARS-H covariance coordinates use the selected pattern. Existing Morton sorting is retained.

The sampling rules are inspired by Section 2 of *The Third Competition on Spatial Statistics for Large Datasets*, DOI 10.1007/s13253-023-00584-9. They are not exact reconstructions of the competition scripts:

- `regular`: independent uniform x,y in [0,1).
- `chessboard`: uniform sampling on even-parity cells of an 8-by-8 board.
- `left-bottom`: square two independent uniform coordinates to concentrate points near zero.
- `clusters`: 85% equal mixture of five Gaussian clusters centered at (0.2,0.2), (0.2,0.8), (0.5,0.5), (0.8,0.2), and (0.8,0.8), plus 15% uniform background. Standard deviations along x are (0.11,0.14,0.16,0.12,0.10), and along y (0.14,0.10,0.12,0.10,0.14). Out-of-domain draws are rejected. Diffuse tails and background replace the earlier isolated disks.

Satellite locations in the paper came from observed data. Load the original coordinates for that case; no synthetic satellite reproduction is claimed. Use a fixed `--seed=42` for repeatability within the same runtime/platform. Standard C rand() is used, so bitwise portability across C libraries is not guaranteed.

Example: append `--location-pattern=chessboard --seed=42` to your usual `Example_Data_Generation` command. The existing generation logger exports locations and observations.

Standalone sampling checks:

```sh
c++ -std=c++17 -I inst/include tests/TestLocationPatternsStandalone.cpp -o test-patterns
./test-patterns
```

These check 100,000 points per pattern, repeatability, bounds, occupancy, concentration, and invalid pattern rejection. Do not compile this assertion-based test with NDEBUG.


## 3D location generation

`LocationPatterns::Fill3D` and the pattern-aware `LocationGenerator::GenerateLocations` overload support 3D regular, chessboard, lower-corner, and diffuse-cluster coordinates. Chessboard cells use parity of all three discretized coordinates. Lower-corner concentration squares independent uniform x,y,z. The 3D cluster mixture contains eight corner components centered at coordinates 0.2 or 0.8 and a center component at (0.5,0.5,0.5), with 15% uniform background. Corner standard deviations are 0.13 + 0.01*(index modulo 3); the center standard deviation is 0.18. Samples outside the unit cube are rejected.

The new 3D sampling checks pass for 90,000 points per pattern. This does not establish PaRSEC covariance support: the public branch's hardware path still explicitly rejects non-default 3D patterns until its backend integration is validated. Do not use these previews as 3D performance measurements. The 2D and 3D extensions are synthetic designs, not exact competition reconstructions.
