# Cable-Polygon Intersection Checker

## Overview

This project processes KMZ files containing cable paths and polygons to determine which cables intersect with specific polygons. It utilizes Shapely, PyKML, and Rtree for spatial analysis.

## Prerequisites

Ensure you have the following installed:

- Python 3.8+
- Required dependencies:
  ```sh
  pip install shapely rtree geopy tqdm pyproj pykml
  ```

## Installation

1. Clone this repository:
   ```sh
   git clone https://github.com/yourusername/cable-intersection-checker.git
   cd cable-intersection-checker
   ```
2. Place your KMZ files (`Patch.kmz` and `polygon_form_placemark_10m.kmz`) in the project directory.

## Usage

Run the script:

```sh
python script.py
```

## Expected Output

- A CSV file (`cable_polygon_intersections.csv`) containing the cable paths, the number of polygons they intersect, and the total cable length in meters.

## Troubleshooting

- If Rtree fails to install, try:
  ```sh
  pip install Rtree
  ```
  Or manually install dependencies using:
  ```sh
  conda install -c conda-forge rtree
  ```
- Ensure your KMZ files contain the correct structure with LineString paths and Polygon elements.

## Contributing

Feel free to fork and contribute! Open an issue if you encounter any problems.
