# GitHub setup quickstart

One-time setup to push this repo to GitHub.

## 1. Unzip and inspect

```bash
unzip green-borders.zip
cd green-borders
ls -la
```

You should see `CLAUDE.md`, `README.md`, `R/`, `docs/`, etc.

## 2. Create GitHub repo (private recommended for thesis)

Option A: via gh CLI

```bash
gh repo create green-borders --private --description "IL cannabis legalization, cross-border entry and crime spillovers"
```

Option B: via web UI at https://github.com/new
- Name: `green-borders`
- Visibility: Private
- Do NOT initialize with README/.gitignore/license (we have our own)

## 3. Initialize and push

```bash
cd green-borders
git init
git add .
git status                    # eyeball that no data/raw/ files are staged
git commit -m "Initial commit: pipeline skeleton and design docs"
git branch -M main
git remote add origin git@github.com:YOUR_USERNAME/green-borders.git
git push -u origin main
```

## 4. Set up local environment

```bash
# Copy Renviron template and fill in your Census API key
cp .Renviron.template .Renviron
# Edit .Renviron to add your CENSUS_API_KEY
```

In RStudio: open `green-borders.Rproj` (or just open the folder), then:

```r
# First-time setup
install.packages("renv")
renv::init(bare = TRUE)
renv::restore()    # installs the pinned package set from renv.lock
```

## 5. Stage raw data (locally, not pushed)

```bash
mkdir -p data/raw/{nibrs,leaic,iowadot,shapefiles}
# Copy your NIBRS Kaplan V9 .rds files into data/raw/nibrs/
# Copy LEAIC .rda into data/raw/leaic/
# Copy Iowa DOT PDFs into data/raw/iowadot/
```

Check `.gitignore` is doing its job:

```bash
git status              # should be clean - no data/raw/ entries
git check-ignore -v data/raw/nibrs/nibrs_offense_segment_2020.rds
# should print: .gitignore:2:data/raw/   data/raw/nibrs/...
```

## 6. First Claude Code task

In your terminal in the repo:

```bash
claude
```

Then issue a tier-0 task from `docs/claude_code_task_examples.md`, e.g.:

> Read CLAUDE.md. Verify the project skeleton matches the directory layout described there. Report any missing files or directories. Do not create files yet.

## 7. Useful housekeeping

```bash
# Add a thesis-defense branch
git checkout -b defense-2026-05-20

# Tag the v1 pipeline snapshot
git tag v1.0-pipeline-skeleton
git push --tags
```
