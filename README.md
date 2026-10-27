# State land-use crop tool

## Download from GitHub

Open a terminal where you want the tool and run:

```text
git clone https://github.com/Yashaskumar123/state_aware_landuse_raster_generator.git
cd state_aware_landuse_raster_generator
```

## What is in `assets`?

- `state_shape.*` contains the state boundaries. Keep all five shapefile files together.
- `plants_lookup.txt` lists the plant names and SWAT+ crop codes you can use.

## What do I do?

1. Put your **original land-use raster** in `user_defined` and name it `landuse.tif`. Another `.tif` or `.tiff` filename also works, but keep only one raster in this folder.
2. Edit `user_defined/state_landuse_crop_lookup.txt` to choose a crop code for each land-use class in each state.
3. Open a terminal in the `state_aware_landuse_raster_generator` folder and run `python reclassify.py`.
4. Find the changed raster in `output/modified_landuse.tif`. The same folder also contains `landuse_change_lookup.csv` (what changed) and `swatplus_landuse_lookup.csv` (raster value to SWAT+ code).

If Python reports a missing package, run `python -m pip install -r requirements.txt` first. For another run, use a new output folder: `python reclassify.py --output-dir output/run2`.

## How do I write the lookup?

The included lookup is ready to edit. Keep **one header** and **one row per state**. Use tabs between entries. After the state ID and index, put each land-use class immediately before its crop code:

```text
State_id	Index	Landuse	change_landuse
MH	1	Kharif	SOYB	Rabi	GRSG	Zaid	MUNG	Double/Triple_Crop	SOYB	Plantation	ORCD	Shifting_Cultivation	PMIL
```

`MH` is Maharashtra; `1` is its state index. `Kharif SOYB` means Kharif cells in Maharashtra become soybean. Choose crop codes from `assets/plants_lookup.txt`. Use the original raster's land-use classes: Kharif=2, Rabi=3, Zaid=4, Double/Triple_Crop=5, Plantation=7, Shifting_Cultivation=13. The output raster value is **original class * 10 + state index**: Maharashtra Kharif becomes **21**.

For this lookup, use `CH` for Chhattisgarh; the tool matches it to the shapefile's `CG` boundary.
