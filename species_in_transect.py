from collections import defaultdict

from shapely.geometry import LineString
from pyproj import Geod

from mnhn_data_cache_client import DataCacheClient

# Define the coordinates,
point1 = (6.106626521825521, 49.62621432351738)
point2 = (6.132649837657024, 49.602114098955234)

# Create a LineString representing the transect,
line = LineString([point1, point2])

# Define a geodetic projection,
geod = Geod(ellps="WGS84")

# Calculate the azimuth of the line,
_, azimuth, _ = geod.inv(point1[0], point1[1], point2[0], point2[1])

# Calculate the midpoint of the line,
midpoint = line.interpolate(0.5, normalized=True)

# Function to create a square polygon rotated by 45° to align with the transect,
def create_rotated_square(center, side_length, azimuth):
    half_diagonal = (side_length / 2) * (2 ** 0.5)  # Diagonal length of the square divided by 2,
    # Adjust the azimuth to account for rotation,
    rotated_azimuths = [azimuth + 45 + i * 90 for i in range(4)]
    # Calculate the corners of the rotated square,
    corners = [geod.fwd(center[0], center[1], a, half_diagonal) for a in rotated_azimuths]
    # Convert to (longitude, latitude) tuples and close the polygon,
    return [[p[0], p[1]] for p in corners] + [[corners[0][0], corners[0][1]]]

# Generate squares starting from the midpoint,
side_length = 200  # in meters,
line_length = geod.line_length(line.xy[0], line.xy[1])
half_squares = int(line_length / (2 * side_length)) + 1  # Number of squares on each side of the midpoint,

squares = []
names = []

for i in range(half_squares):
    # Move outward from the midpoint,
    point_forward = geod.fwd(midpoint.x, midpoint.y, azimuth, i * side_length)
    squares.append(create_rotated_square((point_forward[0], point_forward[1]), side_length, azimuth))
    names.append(f"Transect1-{half_squares+i}")

    if i != 0:    
        point_backward = geod.fwd(midpoint.x, midpoint.y, azimuth + 180, i * side_length)
        squares.insert(0, create_rotated_square((point_backward[0], point_backward[1]), side_length, azimuth))
        names.insert(0, f"Transect1-{half_squares -i}")


dcc = DataCacheClient("serv-data.vm.mnhn.etat.lu")

# Collect statistics per square
sq_stats = []
for sq, sq_name in zip(squares, names):
    sq_stat = defaultdict(lambda: {"vernacular_name": None, "count": 0})
    for query in dcc.search_observations(date_range=["2014-01-01", "2024-12-10"],
                                       polygon=sq):
        hit = query["hits"]["hits"]
        for obs in hit:
            src = obs["_source"]
            sq_stat[src["Taxon_Name"]]["vernacular_name"] = src["Taxon_Common_Names"]
            sq_stat[src["Taxon_Name"]]["count"] += 1
        sq_stat = dict(sorted(sq_stat.items(), key=lambda item: item[1]["count"], reverse=True))
        sq_stats.append(sq_stat)

# Outpout .csv file
print("square_number, count, common_name, species_name, square_link")
for i, sq_stat in enumerate(sq_stats):
    for taxon_name, v in sq_stat.items():
        c = v['count']
        vn = v.get('vernacular_name', 'n/a')
        square_link = 'https://geojson.io/#data=data:application/json,{"type":"FeatureCollection","features":[{"type":"Feature","geometry":{"type":"Polygon","coordinates":[' + repr(squares[i]) + ']}}]}'
        print(f"{i}", end=", ")
        print(f"{c}", end=", ")
        print(f"'{vn}'", end=", ")
        print(f"'{taxon_name}'", end=", ")
        print(f"'{square_link.replace(" ", "")}'", end="")
        print()



