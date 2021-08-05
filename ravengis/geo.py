"""
Tools for performing geospatial translations and transformations.
"""

import json
import logging
from pathlib import Path
from typing import Union

import fiona
from pyproj import CRS
from shapely.geometry import (
    GeometryCollection,
    mapping,
    shape,
)
from shapely.ops import transform

RASTERIO_TIFF_COMPRESSION = "lzw"
LOGGER = logging.getLogger("RavenPy")
WGS84 = 4326


def geom_transform(
    geom: Union[GeometryCollection, shape],
    source_crs: Union[str, int, CRS] = WGS84,
    target_crs: Union[str, int, CRS] = None,
) -> GeometryCollection:
    """Change the projection of a geometry.

    Assuming a geometry's coordinates are in a `source_crs`, compute the new coordinates under the `target_crs`.

    Parameters
    ----------
    geom : Union[GeometryCollection, shape]
      Source geometry.
    source_crs : Union[str, int, CRS]
      Projection identifier (proj4) for the source geometry, e.g. '+proj=longlat +datum=WGS84 +no_defs'.
    target_crs : Union[str, int, CRS]
      Projection identifier (proj4) for the target geometry.

    Returns
    -------
    GeometryCollection
      Reprojected geometry.
    """
    try:
        from functools import partial

        from pyproj import Transformer  # noqa

        if isinstance(source_crs, int or str):
            source = CRS.from_epsg(source_crs)
        else:
            source = source_crs

        if isinstance(target_crs, int or str):
            target = CRS.from_epsg(target_crs)
        else:
            target = target_crs

        transform_func = Transformer.from_crs(source, target, always_xy=True)
        reprojected = transform(transform_func.transform, geom)

        return reprojected
    except Exception as err:
        msg = f"{err}: Failed to reproject geometry"
        LOGGER.error(msg)
        raise Exception(msg)


def generic_vector_reproject(
    vector: Union[str, Path],
    projected: Union[str, Path],
    source_crs: Union[str, CRS] = WGS84,
    target_crs: Union[str, CRS] = None,
) -> None:
    """Reproject all features and layers within a vector file and return a GeoJSON

    Parameters
    ----------
    vector : Union[str, Path]
      Path to a file containing a valid vector layer.
    projected: Union[str, Path]
      Path to a file to be written.
    source_crs : Union[str, dict, CRS]
      Projection identifier (proj4) for the source geometry, Default: '+proj=longlat +datum=WGS84 +no_defs'.
    target_crs : Union[str, dict, CRS]
      Projection identifier (proj4) for the target geometry.

    Returns
    -------
    None
    """

    if target_crs is None:
        raise ValueError("No target CRS is defined.")

    output = {"type": "FeatureCollection", "features": list()}

    if isinstance(vector, Path):
        vector = vector.as_posix()

    for i, layer_name in enumerate(fiona.listlayers(vector)):
        with fiona.open(vector, "r", layer=i) as src:
            with open(projected, "w") as sink:
                for feature in src:
                    # Perform vector reprojection using Shapely on each feature
                    try:
                        geom = shape(feature["geometry"])
                        transformed = geom_transform(geom, source_crs, target_crs)
                        feature["geometry"] = mapping(transformed)
                        output["features"].append(feature)
                    except Exception as err:
                        LOGGER.exception(
                            "{}: Unable to reproject feature {}".format(err, feature)
                        )
                        raise

                sink.write(f"{json.dumps(output)}")
