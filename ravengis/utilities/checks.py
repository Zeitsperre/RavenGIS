"""
Checks for various geospatial and IO conditions.
"""

import logging
import warnings
from collections import abc
from pathlib import Path
from typing import List, Sequence, Tuple, Union

import fiona
import rasterio
from pyproj import CRS
from pyproj.exceptions import CRSError
from shapely.geometry import GeometryCollection, MultiPolygon, Point, shape

from ravengis import io

LOGGER = logging.getLogger("RavenPy")


def boundary_check(
    *args: Sequence[Union[str, Path]],
    max_y: Union[int, float] = 60,
    min_y: Union[int, float] = -60,
) -> None:
    """Verify that boundaries do not exceed specific latitudes for geographic coordinate data. Raise a warning if so.

    Parameters
    ----------
    *args : Sequence[Union[str, Path]]
      listing of strings or paths to files
    max_y : Union[int, float]
      Maximum value allowed for latitude. Default: 60.
    min_y : Union[int, float]
      Minimum value allowed for latitude. Default: -60.
    """
    vectors = (".gml", ".shp", ".geojson", ".gpkg", ".json")
    rasters = (".tif", ".tiff")

    if len(args) == 1 and not isinstance(args[0], str):
        args = args[0]

    for file in args:
        try:
            if str(file).lower().endswith(vectors):
                src = fiona.open(file, "r")
            elif str(file).lower().endswith(rasters):
                src = rasterio.open(file, "r")
            else:
                raise FileNotFoundError()

            try:
                geographic = CRS(src.crs).is_geographic
            except CRSError:
                geographic = True
            src_min_y, src_max_y = src.bounds[1], src.bounds[3]
            if geographic and (src_max_y > max_y or src_min_y < min_y):
                msg = (
                    f"Vector {file} contains geometries in high latitudes."
                    " Verify choice of projected CRS is appropriate for analysis."
                )
                LOGGER.warning(msg)
                warnings.warn(msg, UserWarning)
            if not geographic:
                msg = f"Vector {file} is not in a geographic coordinate system."
                LOGGER.warning(msg)
                warnings.warn(msg, UserWarning)
            src.close()

        except FileNotFoundError:
            msg = f"Unable to read boundaries from {file}"
            LOGGER.error(msg)
            raise
    return


def multipolygon_check(geom: GeometryCollection) -> None:
    """Perform a check to verify a geometry is a MultiPolygon

    Parameters
    ----------
    geom : GeometryCollection

    Returns
    -------
    None
    """
    if not isinstance(geom, GeometryCollection):
        try:
            geom = shape(geom)
        except AttributeError:
            LOGGER.error("Unable to load argument as shapely.geometry.shape().")
            raise

    if isinstance(geom, MultiPolygon):
        LOGGER.warning("Shape is a Multipolygon.")
    return


def feature_contains(
    point: Union[Tuple[Union[int, float, str], Union[str, float, int]], Point],
    shp: Union[str, Path, List[Union[str, Path]]],
) -> Union[dict, bool]:
    """Return the first feature containing a location.

    Parameters
    ----------
    point : Union[Tuple[Union[int, float, str], Union[str, float, int]], Point]
      Geographic coordinates of a point (lon, lat) or a shapely Point.
    shp : Union[str, Path, List[str, Path]]
      Path to the file storing the geometries.

    Returns
    -------
    Union[dict, bool]
      The feature found.

    Notes
    -----
    This is really slow. Another approach is to use the `fiona.Collection.filter` method.
    """

    if isinstance(point, abc.Sequence) and not isinstance(point, str):
        for coord in point:
            if isinstance(coord, (int, float)):
                pass
        point = Point(point)
    elif isinstance(point, Point):
        pass
    else:
        raise ValueError(
            f"point should be shapely.Point or tuple of coordinates, got : {point} of type({type(point)})"
        )

    shape_crs = io.crs_sniffer(shp)

    if isinstance(shp, list):
        shp = shp[0]

    for i, layer_name in enumerate(fiona.listlayers(str(shp))):
        with fiona.open(shp, "r", crs=shape_crs, layer=i) as vector:
            for f in vector.filter(bbox=(point.x, point.y, point.x, point.y)):
                return f

    return False
