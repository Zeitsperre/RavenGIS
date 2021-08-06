import tempfile
from pathlib import Path

import fiona
import numpy as np
import pytest
import rasterio
import shapely.geometry as sgeo

from ravengis import analysis, io, raster
from ravengis.utilities.checks import boundary_check, feature_contains
from ravengis.utilities.testdata import get_local_testdata


class TestOperations:

    zipped_file = get_local_testdata("polygons/mars.zip")
    non_zipped_file = get_local_testdata("polygons/mars.geojson")

    def test_circular_mean_aspect(self):
        northern_angles = np.array([330, 30, 15, 345])
        slight_northeast_angles = np.append(northern_angles, [0.000001])
        eastern_angles = np.arange(45, 125, 1.25)
        southwest_angles = np.array([181, 182.25, 183.5, 222])

        assert analysis.circular_mean_aspect(northern_angles) == 360
        np.testing.assert_almost_equal(
            analysis.circular_mean_aspect(slight_northeast_angles), 0, decimal=3
        )
        assert analysis.circular_mean_aspect(eastern_angles) == 84.375
        np.testing.assert_almost_equal(
            analysis.circular_mean_aspect(southwest_angles), 191.88055987
        )

    def test_address_append(self):
        non_existing_tarred_file = "polygons.tar"

        assert "zip://" in io.address_append(self.zipped_file)
        assert "tar://" in io.address_append(non_existing_tarred_file)
        assert not io.address_append(self.non_zipped_file).startswith(
            ("zip://", "tar://")
        )

    def test_archive_sniffer(self, tmp_path):
        probable_shp = io.archive_sniffer(self.zipped_file)
        assert Path(probable_shp[0]).name == "mars.shp"

        probable_shp = io.archive_sniffer(self.zipped_file, working_dir=tmp_path)
        assert Path(probable_shp[0]).name == "mars.shp"

    def test_archive_extract(self, tmp_path):
        assert self.zipped_file.exists()

        files = list()
        with tempfile.TemporaryDirectory(dir=tmp_path) as tdir:
            files.extend(io.generic_extract_archive(self.zipped_file, output_dir=tdir))
            assert len(files) == 5
            for f in files:
                assert Path(f).exists()
        assert not np.any([Path(f).exists() for f in files])

        files = io.generic_extract_archive(self.zipped_file)
        assert np.all([Path(f).exists() for f in files])


class TestFileInfoFuncs:
    zipped_file = get_local_testdata("polygons/mars.zip")
    geojson_file = get_local_testdata("polygons/mars.geojson")
    raster_file = get_local_testdata(
        "nasa/Mars_MGS_MOLA_DEM_georeferenced_region_compressed.tiff"
    )

    non_existing_file = "unreal.zip"

    def test_raster_datatype_sniffer(self):
        datatype = io.raster_datatype_sniffer(self.raster_file)
        assert datatype.lower() == "uint8"

    def test_crs_sniffer(self):
        assert io.crs_sniffer(self.zipped_file) == 4326
        assert set(io.crs_sniffer(self.geojson_file, self.raster_file)) == {4326}

    def test_boundary_check(self):
        # NOTE: does not presently accept zipped files.

        with pytest.warns(None):
            boundary_check([self.geojson_file, self.raster_file], max_y=80)

        with pytest.warns(UserWarning):
            boundary_check([self.geojson_file, self.raster_file], max_y=15)

        with pytest.raises(FileNotFoundError):
            boundary_check([self.non_existing_file])

    @pytest.mark.skip(reason="Not presently testable")
    def test_multipolygon_check(self):
        pass


class TestGdalOgrFunctions:
    geojson_file = get_local_testdata("polygons/mars.geojson")
    raster_file = get_local_testdata(
        "nasa/Mars_MGS_MOLA_DEM_georeferenced_region_compressed.tiff"
    )

    def test_gdal_aspect_not_projected(self, tmp_path):
        aspect_grid = analysis.gdal_aspect_analysis(self.raster_file)
        np.testing.assert_almost_equal(
            analysis.circular_mean_aspect(aspect_grid), 10.9119033
        )

        # test with creation of a temporary file
        aspect_tempfile = tempfile.NamedTemporaryFile(
            prefix="aspect_", suffix=".tiff", delete=False, dir=tmp_path
        ).name
        aspect_grid = analysis.gdal_aspect_analysis(
            self.raster_file, set_output=aspect_tempfile
        )
        np.testing.assert_almost_equal(
            analysis.circular_mean_aspect(aspect_grid), 10.9119033
        )
        assert Path(aspect_tempfile).stat().st_size > 0

    # Slope values are high due to data values using Geographic CRS
    def test_gdal_slope_not_projected(self, tmp_path):
        slope_grid: np.ndarray = analysis.gdal_slope_analysis(self.raster_file)
        np.testing.assert_almost_equal(slope_grid[:].min(), 0.0)
        np.testing.assert_almost_equal(slope_grid.mean(), 64.4365427)
        np.testing.assert_almost_equal(slope_grid[:].max(), 89.71747, 5)

        slope_tempfile = tempfile.NamedTemporaryFile(
            prefix="slope_", suffix=".tiff", delete=False, dir=tmp_path
        ).name
        slope_grid = analysis.gdal_slope_analysis(
            self.raster_file, set_output=slope_tempfile
        )
        np.testing.assert_almost_equal(slope_grid.mean(), 64.4365427)
        assert Path(slope_tempfile).stat().st_size > 0

    # Slope values are high due to data values using Geographic CRS
    def test_dem_properties(self):
        dem_properties = analysis.dem_prop(self.raster_file)
        np.testing.assert_almost_equal(dem_properties["aspect"], 10.911, 3)
        np.testing.assert_almost_equal(dem_properties["elevation"], 79.0341, 4)
        np.testing.assert_almost_equal(dem_properties["slope"], 64.43654, 5)

        with fiona.open(self.geojson_file) as gj:
            feature = next(iter(gj))
            geom = sgeo.shape(feature["geometry"])

        region_dem_properties = analysis.dem_prop(self.raster_file, geom=geom)
        np.testing.assert_almost_equal(region_dem_properties["aspect"], 280.681, 3)
        np.testing.assert_almost_equal(region_dem_properties["elevation"], 145.8899, 4)
        np.testing.assert_almost_equal(region_dem_properties["slope"], 61.26508, 5)


class TestGenericGeoOperations:
    geojson_file = get_local_testdata("polygons/mars.geojson")
    raster_file = get_local_testdata(
        "nasa/Mars_MGS_MOLA_DEM_georeferenced_region_compressed.tiff"
    )

    def test_raster_warp(self, tmp_path):
        # TODO: It would be awesome if this returned a temporary filepath if no file given.
        # TODO: either use `output` or `reprojected/warped` for these functions.
        reproj_file = tempfile.NamedTemporaryFile(
            prefix="reproj_", suffix=".tiff", delete=False, dir=tmp_path
        ).name
        raster.generic_raster_warp(
            self.raster_file, output=reproj_file, target_crs="EPSG:3348"
        )

        # EPSG:3348 is a very general transformation; Some tolerance should be allowed.
        with rasterio.open(reproj_file) as gt:
            assert gt.crs.to_epsg() == 3348
            np.testing.assert_allclose(gt.bounds.left, -2077535, atol=3)
            np.testing.assert_allclose(gt.bounds.right, 15591620, atol=3)
            np.testing.assert_allclose(gt.bounds.bottom, -4167898, atol=3)
            np.testing.assert_allclose(gt.bounds.top, 5817014, atol=3)

            data = gt.read(1)  # read band 1 (red)
            assert data.min() == 0
            assert data.max() == 255
            np.testing.assert_almost_equal(data.mean(), 60.729, 3)

    def test_warped_raster_slope(self, tmp_path):
        reproj_file = tempfile.NamedTemporaryFile(
            prefix="reproj_", suffix=".tiff", delete=False, dir=tmp_path
        ).name
        raster.generic_raster_warp(
            self.raster_file, output=reproj_file, target_crs="EPSG:3348"
        )
        slope_grid = analysis.gdal_slope_analysis(reproj_file)

        np.testing.assert_almost_equal(slope_grid[:].min(), 0.0)
        np.testing.assert_almost_equal(slope_grid.mean(), 0.0034991)
        np.testing.assert_almost_equal(slope_grid[:].max(), 0.3523546)

    def test_warped_raster_aspect(self, tmp_path):
        reproj_file = tempfile.NamedTemporaryFile(
            prefix="reproj_", suffix=".tiff", delete=False, dir=tmp_path
        ).name
        raster.generic_raster_warp(
            self.raster_file, output=reproj_file, target_crs="EPSG:3348"
        )
        aspect_grid = analysis.gdal_aspect_analysis(reproj_file)

        np.testing.assert_almost_equal(
            analysis.circular_mean_aspect(aspect_grid), 7.780, decimal=3
        )

    def test_raster_clip(self, tmp_path):
        with fiona.open(self.geojson_file) as gj:
            feature = next(iter(gj))
            geom = sgeo.shape(feature["geometry"])

        clipped_file = tempfile.NamedTemporaryFile(
            prefix="clipped_", suffix=".tiff", delete=False, dir=tmp_path
        ).name
        raster.generic_raster_clip(self.raster_file, clipped_file, geometry=geom)

        with rasterio.open(clipped_file) as gt:
            assert gt.crs.to_epsg() == 4326

            data = gt.read(1)  # read band 1 (red)
            assert data.min() == 0
            assert data.max() == 255
            np.testing.assert_almost_equal(data.mean(), 102.8222965)


class TestGIS:
    vector_file = get_local_testdata("polygons/mars.geojson")

    def test_get_bbox_single(self):
        w, s, n, e = io.get_bbox(self.vector_file, all_features=False)
        np.testing.assert_almost_equal(w, -139.8514262)
        np.testing.assert_almost_equal(s, 8.3754794)
        np.testing.assert_almost_equal(n, -117.4753973)
        np.testing.assert_almost_equal(e, 29.6327068)

    def test_get_bbox_all(self):
        w, s, n, e = io.get_bbox(self.vector_file)
        np.testing.assert_almost_equal(w, -139.8514262)
        np.testing.assert_almost_equal(s, 8.3754794)
        np.testing.assert_almost_equal(n, -38.7397456)
        np.testing.assert_almost_equal(e, 64.1757015)

    def test_feature_contains(self):
        point = -69.0, 45
        assert isinstance(feature_contains(point, self.vector_file), dict)
        assert isinstance(
            feature_contains(sgeo.Point(point), self.vector_file), dict
        )
