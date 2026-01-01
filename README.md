# Global warming × nutrient enrichment experiment (R project)

Code and data for the manuscript: **Global intertidal community responses to warming and nutrient enrichment**.

## Repository structure
- `Global multiple stressor project.Rproj`: RStudio project file (sets the working context).
- `Key analyses/`: core manuscript analysis scripts.
- `Additional analyses/`: diagnostics, sensitivity checks, and supplementary processing scripts.
- `Dataframes/`: inputs used by the scripts.

## Quick start (recommended workflow)
1. Install **Git**.
2. Install **Git LFS** (required for the land-cover `.tif`; see below).
3. Clone this repository.
4. Open `Global multiple stressor project.Rproj` in RStudio.
5. Install R packages as needed (each script lists required packages at the top).
6. Run scripts from within the project.

> Note: scripts use the `{here}` package for portable paths (e.g. `here("Dataframes", ...)`).

## Requirements
- R (analyses run on **R 4.4.3**).
- RStudio (recommended).
- Git.
- Git LFS (for land-cover raster; see below).

---

## Large file in this repo (Git LFS) + land-cover input

This repository includes a land-cover raster stored via **Git Large File Storage (Git LFS)**. Git LFS stores a pointer in Git history and fetches the real binary file separately when cloned/pulled.

### How to proceed
1. Install Git LFS (once): https://docs.github.com/en/repositories/working-with-files/managing-large-files/installing-git-large-file-storage
2. In a terminal (or Git Bash), run once:
   - `git lfs install`
3. Clone the repository:
   - `git clone https://github.com/RameshWilson/Global-paper.git`
4. Ensure LFS files are present:
   - `git lfs pull`
5. Open the `.Rproj` in RStudio and run scripts.

### The land-cover file included here
The file used by `Additional analyses/Land use processing_Additional.R` is:

`Dataframes/Land use mapping/lccs_class_2022_v2.1.1_crop_sites_5kmplus500m.tif`

This is a **derived, cropped** GeoTIFF containing only the **`lccs_class`** layer required for the site-buffer regions (crop based on the study extent + margin), produced to reduce file size.

Original source (global product):
- Copernicus Climate Data Store dataset page (includig licence/terms): https://cds.climate.copernicus.eu/datasets/satellite-land-cover

---

## Weather covariates (Visual Crossing)

Weather covariates used in the manuscript were derived from **Visual Crossing** daily weather downloads (retrieved during a time-limited paid subscription).

**Public repository policy:**
- This repository **does not redistribute** any Visual Crossing raw downloads (daily CSVs) (as per VC terms and conditions: https://www.visualcrossing.com/weather-services-terms/).
- This repository also **does not rely on re-downloading** Visual Crossing data to run the core ecological models.
- The script `Additional analyses/Weather data processing_Additional.R` is included **for transparency** (to show how seasonal summaries were computed and subsequently appended into the main dataframe), but it **cannot be run from this public repository alone** as the raw Visual Crossing CSV inputs are not distributed here.

### To reproduce the full weather summaries if desired:
1. Obtain your own weather time series for the same locations and periods (e.g., via Visual Crossing under your own personal use licence, or an alternative provider).
2. Place raw daily files locally (not committed) and run:
   - `Additional analyses/Weather data processing_Additional.R` (with appropriate renaming as required)
3. Merge the resulting per-site seasonal summaries into analysis dataset as needed.

---

## Reproducibility notes
- Paths are relative to the project root (via `{here}`), so analyses should run after cloning.
- Figures are produced interactively in R.

## How to cite
If you use this repository, please cite:
- **Manuscript:** citation + DOI (to be added)
- **Data and code:** hosted by Zenodo: doi.org/10.5281/zenodo.18121701

## Licence
- Code licence: see `LICENSE` (MIT).
- Repository data licence: see `DATA_LICENSE.md` (CC BY 4.0 for original project datasets and author-generated derived tables), **excluding** third-party governed inputs noted below.
- Third-party land-cover data: governed by the Copernicus/ESA licence linked above (not covered by this repo’s licences).
- Weather (Visual Crossing): raw downloads are not redistributed here; use is governed by Visual Crossing terms (not covered by this repo’s licences).

## Contact
- Lead author: Ramesh Wilson
- Email: rameshwilson14@gmail.com
- GitHub: `RameshWilson`
