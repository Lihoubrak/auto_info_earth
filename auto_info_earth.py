import zipfile
from pykml import parser
from shapely.geometry import LineString, Polygon
import os
import csv
from tqdm import tqdm
from pyproj import Geod
from rtree import index  # Import R-tree

# Path to KMZ files
script_dir = os.path.dirname(os.path.abspath(__file__))
cable_kmz_file = os.path.join(script_dir, "Patch.kmz")
polygon_kmz_file = os.path.join(script_dir, "polygon_form_placemark_10m.kmz")

target_folder_name = "PNP Cable"

def extract_kml(kmz_file):
    with zipfile.ZipFile(kmz_file, 'r') as kmz:
        for file_name in kmz.namelist():
            if file_name.endswith('.kml'):
                with kmz.open(file_name) as kml_file:
                    return kml_file.read()
    raise FileNotFoundError(f"No KML file found in {kmz_file}")

cable_kml = parser.fromstring(extract_kml(cable_kmz_file))
polygon_kml = parser.fromstring(extract_kml(polygon_kmz_file))
namespace = "{http://www.opengis.net/kml/2.2}"

# Retrieve the "PNP Cable" folder
kan_cable_folder = next((folder for folder in cable_kml.findall(f".//{namespace}Folder")
                         if folder.find(f"{namespace}name") is not None and folder.find(f"{namespace}name").text.strip() == target_folder_name), None)
if kan_cable_folder is None:
    raise ValueError(f"Folder '{target_folder_name}' not found in KML.")

cable_paths = []
for placemark in kan_cable_folder.findall(f".//{namespace}Placemark"):
    line = placemark.find(f".//{namespace}LineString")
    name = placemark.find(f"{namespace}name").text if placemark.find(f"{namespace}name") is not None else "Unknown"
    if line is not None:
        coords_text = line.find(f".//{namespace}coordinates").text.strip()
        coords = [(float(lon), float(lat)) for lon, lat, *_ in [c.split(',') for c in coords_text.split()]]
        if len(coords) >= 2:
            cable_paths.append((name, LineString(coords)))

# Load polygons and build R-tree index
polygons = []
polygon_names = []
polygon_bboxes = []
rtree_idx = index.Index()

print("\n🔄 Loading polygons into R-tree...")
for i, placemark in enumerate(tqdm(polygon_kml.findall(f".//{namespace}Placemark"), desc="📍 Processing polygons", unit="polygon")):
    polygon = placemark.find(f".//{namespace}Polygon")
    name = placemark.find(f"{namespace}name").text if placemark.find(f"{namespace}name") is not None else "Polygon"
    if polygon is not None:
        coords_text = polygon.find(f".//{namespace}outerBoundaryIs//{namespace}LinearRing//{namespace}coordinates").text.strip()
        coords = [(float(lon), float(lat)) for lon, lat, *_ in [c.split(',') for c in coords_text.split()]]
        if len(coords) >= 3:
            poly = Polygon(coords)
            bbox = poly.bounds  # (minx, miny, maxx, maxy)
            polygons.append(poly)
            polygon_names.append(name)
            polygon_bboxes.append(bbox)
            rtree_idx.insert(i, bbox)  # Add to R-tree

def calculate_real_length(line):
    geod = Geod(ellps="WGS84")
    total_length = 0
    coords = list(line.coords)
    
    for i in range(len(coords) - 1):
        lon1, lat1 = coords[i]
        lon2, lat2 = coords[i + 1]
        _, _, distance = geod.inv(lon1, lat1, lon2, lat2)
        total_length += distance

    return total_length  # Return total length in meters

# Lists to store results
filtered_results = []  # For cables that intersect with polygons
non_intersecting_results = []  # For cables that do not intersect with any polygons

print("\n🔄 Checking cable paths intersecting with polygons...")

for name, path in tqdm(cable_paths, desc="🔍 Processing cables", unit="cable"):
    path_bbox = path.bounds  # (minx, miny, maxx, maxy)
    
    # Query R-tree index to find polygons whose bounding box intersects with the cable path
    candidate_indices = list(rtree_idx.intersection(path_bbox))
    
    crossed_polygons = [polygon_names[i] for i in candidate_indices if path.crosses(polygons[i]) or path.intersects(polygons[i])]
    
    real_length_meters = calculate_real_length(path)
    
    if crossed_polygons:
        num_crossed = len(crossed_polygons)
        filtered_results.append([name, num_crossed, ", ".join(crossed_polygons), real_length_meters])
    else:
        non_intersecting_results.append([name, real_length_meters])

# Write results for intersecting cables
csv_file_intersecting = os.path.join(script_dir, "cable_polygon_intersections.csv")
with open(csv_file_intersecting, mode='w', newline='') as file:
    writer = csv.writer(file)
    writer.writerow(["Cable Path", "Polygons Crossed", "Polygon Names", "Total Length (m)"])
    writer.writerows(filtered_results)

# Write results for non-intersecting cables
csv_file_non_intersecting = os.path.join(script_dir, "cable_no_polygon_intersections.csv")
with open(csv_file_non_intersecting, mode='w', newline='') as file:
    writer = csv.writer(file)
    writer.writerow(["Cable Path", "Total Length (m)"])
    writer.writerows(non_intersecting_results)

print(f"\n✅ Completed!")
print(f"Results for intersecting cables saved at: {csv_file_intersecting}")
print(f"Results for non-intersecting cables saved at: {csv_file_non_intersecting}")