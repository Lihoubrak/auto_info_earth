import zipfile
from pykml import parser
from shapely.geometry import LineString, Polygon
import os
import csv
from tqdm import tqdm
from pyproj import Geod
from rtree import index  # Import R-tree

# Đường dẫn đến file KMZ
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

# Lấy thư mục "SPE Cable"
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

# Load polygons và xây dựng chỉ mục R-tree
polygons = []
polygon_names = []
polygon_bboxes = []
rtree_idx = index.Index()

print("\n🔄 Đang tải các đa giác vào R-tree...")
for i, placemark in enumerate(tqdm(polygon_kml.findall(f".//{namespace}Placemark"), desc="📍 Đang xử lý đa giác", unit="polygon")):
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
            rtree_idx.insert(i, bbox)  # Thêm vào R-tree

def calculate_real_length(line):
    geod = Geod(ellps="WGS84")
    total_length = 0
    coords = list(line.coords)
    
    for i in range(len(coords) - 1):
        lon1, lat1 = coords[i]
        lon2, lat2 = coords[i + 1]
        _, _, distance = geod.inv(lon1, lat1, lon2, lat2)
        total_length += distance

    return total_length  # Trả về tổng chiều dài theo mét

filtered_results = []
print("\n🔄 Đang kiểm tra các đường cáp giao cắt với đa giác...")

for name, path in tqdm(cable_paths, desc="🔍 Đang xử lý cáp", unit="cable"):
    path_bbox = path.bounds  # (minx, miny, maxx, maxy)
    
    # Truy vấn chỉ mục R-tree để tìm các polygon có bounding box giao nhau với đường cáp
    candidate_indices = list(rtree_idx.intersection(path_bbox))
    
    crossed_polygons = [polygon_names[i] for i in candidate_indices if path.crosses(polygons[i]) or path.intersects(polygons[i])]
    
    if crossed_polygons:
        num_crossed = len(crossed_polygons)
        real_length_meters = calculate_real_length(path)
        filtered_results.append([name, num_crossed, ", ".join(crossed_polygons), real_length_meters])

csv_file = os.path.join(script_dir, "cable_polygon_intersections.csv")
with open(csv_file, mode='w', newline='') as file:
    writer = csv.writer(file)
    writer.writerow(["Cable Path", "Polygons Crossed", "Polygon Names", "Total Length (m)"])
    writer.writerows(filtered_results)

print(f"\n✅ Hoàn thành! Kết quả đã lưu tại: {csv_file}")
