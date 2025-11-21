# 1.0.20251001
- Beta release of OCN package.

# 1.1.20251005
- added support for periodic boundary conditions
- changed how intermediate grid states are saved and returned during optimization
- updated the API for requesting intemediate grid states during optimization
- added support for custom cooling schedules

# 1.2.20251010
- api changes for better usability
- added unit testing
- fixed history tracking bug
- fixed a bug in the rng algorithm that impacted reproducibility in some cases
- odd-sized grids are now supported
- fixed a bug preventing early exit in some situations during optimization when the convergence criterion is met

# 1.3.20251010
- fixed critical bug preventing export to xarray datasets

# 1.3.20251011
- fixed a bug where optimization did not work for single-outlet networks

# 1.4.20251029
- Added convenience functions for parallel execution using multiprocessing
- Minor bugfixes

# 1.5.20251121
- Added elevation calculation
- Removed watershed ID attribute from exported rasters (now includes elevation instead). You can still compute watershed IDs on the graph representation if needed.
- Updated demo notebook
- Minor bugfixes