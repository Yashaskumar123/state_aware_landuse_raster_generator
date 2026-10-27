"""Apply state-specific crop codes to a categorical land-use GeoTIFF."""

from __future__ import annotations

import argparse
import csv
import math
from collections import Counter
from pathlib import Path

import numpy as np
import rasterio
import shapefile
from rasterio.crs import CRS
from rasterio.enums import Resampling
from rasterio.features import bounds, geometry_mask
from rasterio.warp import transform_geom
from rasterio.windows import Window, from_bounds


ROOT = Path(__file__).resolve().parent
CLASS_NAMES = {
    2: "Kharif",
    3: "Rabi",
    4: "Zaid",
    5: "Double/Triple_Crop",
    7: "Plantation",
    13: "Shifting_Cultivation",
}
NAME_TO_VALUE = {
    "kharif": 2,
    "kharifcrop": 2,
    "rabi": 3,
    "rabicrop": 3,
    "zaid": 4,
    "zaidcrop": 4,
    "doubletriplecrop": 5,
    "doublecrop": 5,
    "triplecrop": 5,
    "plantation": 7,
    "shiftingcultivation": 13,
}
# These unchanged codes come from the original Godavari land-use DBF mapping.
# WATR and URBN are SWAT+ special codes; the other entries occur in plants.plt.
UNCHANGED_CODES = {
    1: "URBN",
    6: "BSVG",
    8: "FRSE",
    9: "FRSD",
    10: "RNGB",
    12: "RNGE",
    14: "BSVG",
    16: "WATR",
    17: "WATR",
}
CHANGE_COLUMNS = (
    "LANDUSE_ID", "SWAT_CODE", "STATE_ID", "STATE_INDEX", "ORIGINAL_VALUE",
    "ORIGINAL_CLASS", "PLANT", "PIXEL_COUNT",
)


def class_value(label: str) -> int:
    clean = "".join(char for char in label.lower() if char.isalnum())
    if clean.isdecimal() and int(clean) in CLASS_NAMES:
        return int(clean)
    if clean not in NAME_TO_VALUE:
        raise ValueError(f"Unknown land-use class {label!r}; use 2, 3, 4, 5, 7, or 13")
    return NAME_TO_VALUE[clean]


def read_plants(path: Path) -> dict[str, str]:
    plants = {}
    with path.open("r", encoding="utf-8-sig", newline="") as stream:
        reader = csv.DictReader(stream, delimiter="\t")
        if reader.fieldnames != ["ID", "plant", "code"]:
            raise ValueError(f"{path}: expected tab-separated ID, plant, code columns")
        for row in reader:
            code = row["code"].strip().upper()
            if not code or not row["plant"].strip():
                raise ValueError(f"{path}: incomplete plant record")
            plants[code] = row["plant"].strip()
    return plants


def read_changes(path: Path, plants: dict[str, str]) -> list[dict]:
    changes = []
    used_states, used_indexes, used_pairs = set(), set(), set()
    with path.open("r", encoding="utf-8-sig", newline="") as stream:
        reader = csv.reader(stream, delimiter="\t")
        header = next(reader, None)
        if header != ["State_id", "Index", "Landuse", "change_landuse"]:
            raise ValueError(f"{path}: header must be State_id, Index, Landuse, change_landuse")
        for line_number, cells in enumerate(reader, 2):
            if not cells or all(not cell.strip() for cell in cells):
                continue
            if len(cells) < 4 or (len(cells) - 2) % 2:
                raise ValueError(f"{path}:{line_number}: expected state, index, then land-use/code pairs")
            state = cells[0].strip().upper()
            # The user's CH means Chhattisgarh; state_shape.shp calls it CG.
            shape_state = "CG" if state == "CH" else state
            if not state:
                raise ValueError(f"{path}:{line_number}: state ID is empty")
            try:
                index = int(cells[1].strip())
            except ValueError as exc:
                raise ValueError(f"{path}:{line_number}: index must be an integer") from exc
            if not 1 <= index <= 9:
                raise ValueError(f"{path}:{line_number}: index must be 1 through 9 for value = class*10+index")
            if shape_state in used_states:
                raise ValueError(f"{path}:{line_number}: repeated state {state}")
            if index in used_indexes:
                raise ValueError(f"{path}:{line_number}: repeated index {index}")
            used_states.add(shape_state)
            used_indexes.add(index)
            for position in range(2, len(cells), 2):
                original = class_value(cells[position].strip())
                code = cells[position + 1].strip().upper()
                if code not in plants:
                    raise ValueError(f"{path}:{line_number}: {code!r} is absent from plants_lookup.txt")
                if (shape_state, original) in used_pairs:
                    raise ValueError(f"{path}:{line_number}: repeated class {original} for {state}")
                used_pairs.add((shape_state, original))
                changes.append({
                    "state": state,
                    "shape_state": shape_state,
                    "index": index,
                    "original": original,
                    "new": original * 10 + index,
                    "code": code,
                    "plant": plants[code],
                })
    if not changes:
        raise ValueError(f"{path}: no changes were defined")
    return changes


def read_geometries(shp: Path, raster_crs: CRS, required_states: set[str]) -> dict:
    prj = shp.with_suffix(".prj")
    if not prj.exists():
        raise ValueError(f"Missing state projection file: {prj}")
    shape_crs = CRS.from_wkt(prj.read_text(encoding="utf-8"))
    geometries = {}
    with shapefile.Reader(str(shp)) as reader:
        for shape_record in reader.iterShapeRecords():
            state = str(shape_record.record["State"]).strip().upper()
            if state not in required_states:
                continue
            if state in geometries:
                raise ValueError(f"State {state} occurs multiple times in {shp}")
            geom = shape_record.shape.__geo_interface__
            if shape_crs != raster_crs:
                geom = transform_geom(shape_crs, raster_crs, geom)
            geometries[state] = geom
    missing = required_states - geometries.keys()
    if missing:
        raise ValueError(f"State(s) missing from {shp}: {', '.join(sorted(missing))}")
    return geometries


def clipped_window(geom: dict, dataset) -> Window | None:
    candidate = from_bounds(*bounds(geom), transform=dataset.transform)
    col_start = max(0, math.floor(candidate.col_off))
    row_start = max(0, math.floor(candidate.row_off))
    col_stop = min(dataset.width, math.ceil(candidate.col_off + candidate.width))
    row_stop = min(dataset.height, math.ceil(candidate.row_off + candidate.height))
    if col_start >= col_stop or row_start >= row_stop:
        return None
    return Window(col_start, row_start, col_stop - col_start, row_stop - row_start)


def tiles(window: Window, size: int = 1024):
    row_stop = int(window.row_off + window.height)
    col_stop = int(window.col_off + window.width)
    for row in range(int(window.row_off), row_stop, size):
        for col in range(int(window.col_off), col_stop, size):
            yield Window(col, row, min(size, col_stop - col), min(size, row_stop - row))


def tally(array: np.ndarray, counter: Counter) -> None:
    values, counts = np.unique(array[array > 0], return_counts=True)
    counter.update({int(value): int(count) for value, count in zip(values, counts)})


def process(raster: Path, changes: list[dict], states: Path, change_lookup: Path,
            output_dir: Path) -> tuple[Path, Path, Path]:
    output_dir.mkdir(parents=True, exist_ok=True)
    output_raster = output_dir / "modified_landuse.tif"
    change_csv = output_dir / "landuse_change_lookup.csv"
    swat_csv = output_dir / "swatplus_landuse_lookup.csv"
    for path in (output_raster, change_csv, swat_csv):
        if path.exists():
            raise FileExistsError(f"Output already exists: {path}; choose a different --output-dir")
    if raster.resolve() == output_raster.resolve():
        raise ValueError("Input and output raster paths must differ")

    changes_by_state = {}
    for change in changes:
        changes_by_state.setdefault(change["shape_state"], []).append(change)
    original_counts, assigned_counts, output_counts = Counter(), Counter(), Counter()
    with rasterio.open(raster) as src:
        if src.crs is None:
            raise ValueError(f"Input raster has no coordinate system: {raster}")
        if src.count != 1 or not np.issubdtype(np.dtype(src.dtypes[0]), np.integer):
            raise ValueError("Input must be a single-band integer categorical raster")
        if src.transform.b or src.transform.d or src.transform.a <= 0 or src.transform.e >= 0:
            raise ValueError("Input raster must be north-up and have no rotation")
        geometries = read_geometries(states, src.crs, set(changes_by_state))
        largest = max(np.iinfo(np.dtype(src.dtypes[0])).max, max(c["new"] for c in changes))
        dtype = "uint16" if largest <= 65535 else "uint32"
        if largest > np.iinfo(np.uint32).max:
            raise ValueError("Land-use IDs are too large for an unsigned 32-bit raster")
        profile = src.profile.copy()
        profile.update(dtype=dtype, nodata=0, count=1, compress="deflate", predictor=2,
                       tiled=True, blockxsize=256, blockysize=256, BIGTIFF="IF_SAFER")
        try:
            print("Copying source raster...", flush=True)
            with rasterio.open(output_raster, "w", **profile) as dst:
                for _, window in src.block_windows(1):
                    source = src.read(1, window=window)
                    valid = (src.read_masks(1, window=window) != 0) & (source > 0)
                    block = np.where(valid, source, 0).astype(dtype)
                    tally(block, original_counts)
                    dst.write(block, 1, window=window)
                collisions = sorted(set(original_counts) & {c["new"] for c in changes})
                if collisions:
                    raise ValueError(
                        f"New class values already occur in the input: {collisions}. "
                        "Use an unreclassified source raster."
                    )
                dst.update_tags(VALUE_RULE="original_value * 10 + lookup_Index",
                                STATE_LOOKUP=change_lookup.name,
                                AREA_OR_POINT="Area")
                dst.update_tags(1, RepresentationType="THEMATIC")
            print("Applying state crop changes...", flush=True)
            with rasterio.open(output_raster, "r+") as dst:
                for state, state_changes in changes_by_state.items():
                    geometry = geometries[state]
                    extent = clipped_window(geometry, src)
                    if extent is None:
                        continue
                    for window in tiles(extent):
                        source = src.read(1, window=window)
                        if not np.isin(source, [c["original"] for c in state_changes]).any():
                            continue
                        valid = src.read_masks(1, window=window) != 0
                        inside = geometry_mask([geometry], out_shape=source.shape,
                                               transform=src.window_transform(window),
                                               invert=True, all_touched=False)
                        if not inside.any():
                            continue
                        output = dst.read(1, window=window)
                        modified = False
                        for change in state_changes:
                            selected = valid & inside & (source == change["original"])
                            count = int(np.count_nonzero(selected))
                            if count:
                                output[selected] = change["new"]
                                assigned_counts[change["new"]] += count
                                modified = True
                        if modified:
                            dst.write(output, 1, window=window)
                levels = [scale for scale in (2, 4, 8, 16, 32)
                          if src.width // scale >= 1 and src.height // scale >= 1]
                if levels:
                    dst.build_overviews(levels, Resampling.nearest)
                    dst.update_tags(ns="rio_overview", resampling="nearest")
            print("Checking output and writing lookups...", flush=True)
            with rasterio.open(output_raster) as result:
                for _, window in result.block_windows(1):
                    tally(result.read(1, window=window), output_counts)
            if sum(original_counts.values()) != sum(output_counts.values()):
                raise RuntimeError("Positive pixel count changed unexpectedly")
            for change in changes:
                if output_counts[change["new"]] != assigned_counts[change["new"]]:
                    raise RuntimeError(f"Overlapping states or inconsistent output at {change['new']}")
            for value, count in original_counts.items():
                if value not in CLASS_NAMES and output_counts[value] != count:
                    raise RuntimeError(f"Untargeted class {value} changed unexpectedly")
            if not sum(assigned_counts.values()):
                raise RuntimeError("No raster cells matched the state/crop lookup")

            with change_csv.open("w", newline="", encoding="utf-8") as stream:
                writer = csv.DictWriter(stream, fieldnames=CHANGE_COLUMNS)
                writer.writeheader()
                for change in sorted(changes, key=lambda c: (c["index"], c["original"])):
                    writer.writerow({
                        "LANDUSE_ID": change["new"],
                        "SWAT_CODE": change["code"],
                        "STATE_ID": change["state"],
                        "STATE_INDEX": change["index"],
                        "ORIGINAL_VALUE": change["original"],
                        "ORIGINAL_CLASS": CLASS_NAMES[change["original"]],
                        "PLANT": change["plant"],
                        "PIXEL_COUNT": output_counts[change["new"]],
                    })
            assigned_codes = {change["new"]: change["code"] for change in changes}
            with swat_csv.open("w", newline="", encoding="utf-8") as stream:
                writer = csv.writer(stream)
                writer.writerow(["LANDUSE_ID", "SWAT_CODE"])
                for value in sorted(output_counts):
                    writer.writerow([value, assigned_codes.get(value, UNCHANGED_CODES.get(value, ""))])
            unknown = sorted(v for v in output_counts if v not in assigned_codes and v not in UNCHANGED_CODES)
            print(f"Created raster: {output_raster}")
            print(f"Created change lookup: {change_csv}")
            print(f"Reclassified cells: {sum(assigned_counts.values()):,}")
            print(f"Positive cells preserved: {sum(output_counts.values()):,}")
            if unknown:
                print(f"WARNING: No SWAT code is defined for unchanged input value(s) {unknown}; "
                      f"fill their blank codes in {swat_csv} before importing it into SWAT+.")
            return output_raster, change_csv, swat_csv
        except Exception:
            # An interrupted or invalid run must not leave a usable-looking raster.
            for path in (output_raster, change_csv, swat_csv):
                if path.exists():
                    path.unlink()
            raise


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--raster", type=Path, help="Input categorical GeoTIFF; default: only TIFF in user_defined")
    parser.add_argument("--changes", type=Path,
                        default=ROOT / "user_defined" / "state_landuse_crop_lookup.txt")
    parser.add_argument("--states", type=Path, default=ROOT / "assets" / "state_shape.shp")
    parser.add_argument("--plants", type=Path, default=ROOT / "assets" / "plants_lookup.txt")
    parser.add_argument("--output-dir", type=Path, default=ROOT / "output")
    args = parser.parse_args()
    if args.raster is None:
        candidates = sorted((ROOT / "user_defined").glob("*.tif"))
        candidates += sorted((ROOT / "user_defined").glob("*.tiff"))
        if len(candidates) != 1:
            parser.error(f"Expected exactly one TIFF in {ROOT / 'user_defined'}; found {len(candidates)}. "
                         "Add one raster or pass --raster PATH.")
        args.raster = candidates[0]
    try:
        plants = read_plants(args.plants)
        changes = read_changes(args.changes, plants)
        process(args.raster, changes, args.states, args.changes, args.output_dir)
    except (ValueError, FileNotFoundError, FileExistsError, RuntimeError) as exc:
        parser.exit(1, f"Error: {exc}\n")


if __name__ == "__main__":
    main()
