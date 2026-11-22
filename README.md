# Oncoprint Builder

A small R tool that builds **oncoprints** (oncoplots) from annotated variant data, using two plain-text files (a whitelist and a blacklist) and a reviewable CSV to control which genes go in - instead of hardcoding a gene panel inside an R script every time you want to change the figure.

If you've ever built an oncoprint by editing `genes <- c("TP53", "ATRX", ...)` in R, re-running the whole script to see what changed, and then losing that panel when you moved to the next cohort, this replaces that with: edit a text file → run a script → review a CSV → get your figure.

It takes either a long annotated variant table you already have, or raw [SnpEff](https://pcingola.github.io/SnpEff/)-annotated VCFs, and renders with [ComplexHeatmap](https://jokergoo.github.io/ComplexHeatmap-reference/book/)'s `oncoPrint()`.

## Table of contents

- [What this actually does](#what-this-actually-does)
- [Requirements](#requirements)
- [Files in this repo](#files-in-this-repo)
- [Quick start](#quick-start)
- [Full walkthrough](#full-walkthrough)
- [Note file format (`gene_whitelist.txt` / `gene_blacklist.txt`)](#note-file-format-gene_whitelisttxt--gene_blacklisttxt)
- [`consequence_classes.csv` columns](#consequence_classescsv-columns)
- [`gene_list.csv` columns](#gene_listcsv-columns)
- [Input table format](#input-table-format)
- [Configuration reference](#configuration-reference)
- [Output structure](#output-structure)
- [Troubleshooting](#troubleshooting)
- [Citing ComplexHeatmap / SnpEff / vcfR](#citing-complexheatmap--snpeff--vcfr)
- [Citing this repository](#citing-this-repository)
- [License](#license)
- [Provenance](#provenance)

## What this actually does

An oncoprint is three decisions stacked on top of each other, and in most scripts all three are buried in the same block of R:

1. **Which variants count as an alteration at all** - is a synonymous change an alteration? an intron variant? a `MODIFIER`-impact upstream variant?
2. **Which consequence terms collapse into which coloured category** - SnpEff emits ~40 distinct terms; an oncoprint shows maybe five or six colours.
3. **Which genes are rows** - the part you actually argue about, and the part most likely to change between drafts.

This tool pulls each one out into a file you can edit and look at:

| Decision | Lives in |
|---|---|
| What counts as an alteration | `impact_whitelist` in the CONFIG chunk (coarse, changes rarely) |
| Consequence → colour category | `consequence_classes.csv` (one row per SnpEff term) |
| Which genes are rows | `gene_whitelist.txt` + `gene_blacklist.txt`, then `gene_list.csv` |

The gene decision runs in two stages, the same way [GMT Collection Builder](https://github.com/navinvent/GMT_Collection_Builder) handles pathways:

1. **Pattern filter** - a whitelist of gene symbols/patterns to consider and a blacklist to always drop (long genes, HLA, pseudogenes - the usual passenger noise). The blacklist wins.
2. **Recurrence filter** - of the genes that survive, keep those altered in at least *n* samples, up to a cap.

The result is written to `gene_list.csv`, which you review and hand-edit - flip genes off, delete rows, or add a gene by name that no pattern happened to catch - and the figure is built from *that CSV*, not from the filters directly. A one-off manual decision ("I want RB1 in this figure even though it's only in one sample") survives independently of whatever the whitelist and blacklist say, and is still there, with your note next to it, when you come back in three months.

## Requirements

- R (any recent version; developed against 4.6)
- From CRAN:

```r
install.packages(c("dplyr", "tidyr", "stringr", "vcfR"))
```

- From Bioconductor:

```r
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install("ComplexHeatmap")
```

Optional, picked up automatically if present:

- **`data.table`** - much faster reading of large variant tables. With `R.utils` alongside it, that speedup also applies to `.csv.gz`.
- **`qs`** - only needed if your variant table is a `.qs` file.
- **`vcfR`** is only needed for Section 0 (parsing VCFs). Skip it if you already have a variant table.

No internet access is needed at run time.

## Files in this repo

| File | Purpose |
|---|---|
| `Build_Oncoprint.rmd` | The script. Open it in RStudio/Positron and run its chunks in order (see below). |
| `gene_whitelist.txt` | Gene symbols/patterns to consider, one per line. |
| `gene_blacklist.txt` | Gene symbols/patterns to always drop, one per line - wins over the whitelist. |
| `consequence_classes.csv` | Maps each SnpEff consequence term to a colour category, and says which categories are drawn. Edit freely. |
| `gene_list.example.csv` | A **static preview** of what the working CSV looks like. Illustrative - not read by the script, and not real results. |
| `sample_annotation.example.csv` | A **static preview** of the optional per-sample metadata file. Also illustrative. |
| `gene_list.csv` | Generated by Section 1, hand-edited by you, consumed by Section 2. The actual source of truth for which genes are rows. Gitignored - a per-run working file, regenerate it yourself. |
| `sample_annotation.csv` | Optional, written by you, read by Section 2 for the annotation bars. Gitignored - it's your cohort's metadata, not the tool's. |
| `output/` | Where figures and tables land. Gitignored, regenerated on every run. |

## Quick start

1. Clone this repo.
2. Open `Build_Oncoprint.rmd` in RStudio or Positron.
3. In the **CONFIG chunk** near the top, change `base_dir` to wherever you cloned this repo (it's an absolute path - there's no way for the script to guess this). Then point it at your data: either set `variant_table` to an annotated variant table you already have, or set `vcf_dir` to a folder of SnpEff-annotated VCFs.
4. Run the CONFIG chunk.
5. **If you set `vcf_dir`:** run **Section 0**. It parses the VCFs into `output/variants_annotated.csv.gz`, which is where `variant_table` already points.
6. Run **Section 1**. This applies the filters and writes `gene_list.csv`.
7. Open `gene_list.csv` and review it. Edit freely.
8. Run **Section 2**. This reads `gene_list.csv` and writes the oncoprint to `output/`.

That's the whole loop. Steps 7-8 repeat any time you just want to change which genes are rows, without re-running Sections 0 or 1.

## Full walkthrough

### Step 0: get your variants into a long table

Section 2 and Section 1 both work from **one row per variant, per gene, per consequence**. If a variant hits two genes, that's two rows. If SnpEff gave it two consequences, that's two rows.

If you already have such a table - from your own pipeline, from `vcfR`, from anything - skip Section 0 and point `variant_table` at it, then tell the script which of your columns mean what using `cols` in CONFIG. Only `sample`, `gene` and `consequence` are required.

If you have raw SnpEff VCFs, set `vcf_dir` and run Section 0. It:

1. Finds every `.vcf` / `.vcf.gz` under that folder, recursively.
2. Drops variants below `min_qual` (default 30). A variant whose `QUAL` is missing or non-numeric - `.`, which Mutect2 and several other callers write - cannot be compared to a threshold, so `min_qual` is skipped for that file and the skip is reported. It is not treated as a failed filter.
3. Pulls the `ANN=` field out of `INFO` and splits it - SnpEff packs multiple annotations per variant, comma-separated, each one 16 pipe-separated fields.
4. Splits `&`-joined consequence terms (`splice_acceptor_variant&intron_variant`) into separate rows.
5. Derives per-sample metadata from each **file name** using `filename_fields`, which defaults to the `<Sample>_<CancerType>_<State>.<caller>.ann.vcf` convention. Change those regexes to match your own naming.
6. Writes `output/variants_annotated.csv.gz`.

Reading a cohort of exome VCFs is the slow part of this tool - minutes, not seconds, and memory-hungry. It's a separate section precisely so you only do it once.

### Step 1: decide what counts as an alteration

Two places, deliberately:

**`impact_whitelist` in CONFIG** - SnpEff assigns every annotation an impact of `HIGH`, `MODERATE`, `LOW` or `MODIFIER`. The default keeps `HIGH` and `MODERATE`, which is the usual "protein-altering" cut. Set it to `NULL` to keep everything and let `consequence_classes.csv` do all the work.

**`consequence_classes.csv`** - one row per SnpEff term, mapping it to a `class` (the thing that gets a colour in the legend) and a `plot` flag (whether that class is drawn at all). The defaults draw nine classes - Amplification, Deletion, Structural, Fusion, Nonsense, Frameshift, Splice, Inframe_indel and Missense - and carry Synonymous / UTR / Intron / Noncoding / Other as `plot=FALSE` rows so they're one edit away rather than missing. The four CNV and structural classes only appear if your variant table actually uses those consequence terms; SnpEff alone will never emit them.

Any consequence term in your data that has **no row in this file is dropped and reported by name** at the end of Section 1. That report is the point: it tells you exactly what you'd be throwing away, instead of silently binning it.

### Step 1b: filter out inherited variants

**Read this before you trust the figure.** The gene panel is chosen by recurrence across samples, and an inherited variant is recurrent in every sample that carries it by definition. Point this tool at tumour-only VCFs with no germline gate and the top rows can be ordinary population polymorphisms that merely cleared the impact filter. On one cohort where it could be checked directly against a matched normal, **94% of the variants clearing `QUAL >= 30` plus HIGH/MODERATE impact were also present in that donor's own germline sample.**

Dropping the germline *samples* from your cohort does not fix this. That removes the normals; it leaves every inherited variant sitting inside the tumour samples untouched. The normals are the evidence, not the contamination.

Three gates, in section 3b of CONFIG. Each is **off until you configure it**, and Section 1 prints which ones fired, so a misconfigured gate cannot pass for a working one.

**Gate 1 - subtract variants seen in normals.** Set `cols$type` to the column that marks sample type and `germline_type_value` to the value meaning "normal" (`"Germline"` for the default filename convention). Needs `cols$chrom`, `cols$pos`, `cols$ref` and `cols$alt` so a variant is matched exactly; position alone would conflate a hotspot with a benign change at the same base, and the tool refuses rather than doing that.

`germline_scope` is the honesty setting:

| Scope | What it does | When it's right |
|---|---|---|
| `"matched"` | Subtracts a normal only from its own donor's samples. Needs `cols$donor`. | Nearly always. It is all a handful of normals can honestly support. |
| `"cohort"` | Treats every normal as a panel of normals for every sample. | Only with a real panel. One individual is hom-ref at a large share of common SNPs, so a one-donor "panel" both misses most of them and injects that person's private variants into everyone else's blacklist. Section 1 warns below five normals. |

**Gate 2 - population allele frequency.** Set `cols$pop_af` to a population-frequency column and `max_pop_af` to the fraction above which a variant is called inherited (default `0.01`, following GATK `FilterMutectCalls`' `--max-population-af`; there is no consensus value and it does not transfer between assays, so treat it as yours to set).

> **The trap worth reading twice.** freebayes writes its own `INFO` key named `AF`, and it is the **within-sample allele fraction** - about 0.50 for a het, 1.00 for a hom - not population frequency. Gating on it would delete nearly every real variant while reading like a textbook filter. A column named exactly `AF` is **refused**, not trusted. Annotate with an explicit prefix so the two can never be confused:
>
> ```
> SnpSift annotate -info AF,AF_grpmax -name gnomad_ gnomad.exomes.v4.1.sites.vcf.bgz in.vcf > out.vcf
> ```
>
> then set `cols$pop_af <- "gnomad_AF_grpmax"`. The `-name gnomad_` prefix is mandatory, not style: without it SnpSift overwrites freebayes' own `AF` on matching rows only, leaving one column whose meaning varies per row.

A variant with no population-frequency value is **kept**. Absence from a population database is evidence of rarity, and dropping what you cannot score would quietly delete the rare pathogenic variants this gate exists to preserve.

**Gate 3 - a flag, never a filter.** With `cols$ao` and `cols$ro` set (freebayes `AO`/`RO`), `gene_list.csv` gains a `median_vaf` column. A germline het sits near 0.50 and a hom near 1.00 in every carrier, while a somatic variant is usually subclonal and variable. It is a proxy confounded by tumour purity, ploidy and loss of heterozygosity, so it is reported for your review and never used to drop a row.

**Cohort recurrence is deliberately not offered.** Flagging a variant because it recurs across most of the cohort cannot distinguish a polymorphism from the recurrent hotspot that defines a histology. In a seven-sample DMG cohort, H3-3A K27M occurs at the identical base in roughly 75% of cases, so a "present in every sample" rule deletes the single most important row about one time in eight. At realistic cohort sizes the estimator has no resolving power anyway: an observed 50% at n = 14 has a 95% interval of 23% to 77%. Gates 1 and 2 rest on evidence about the variant; cohort frequency rests on how often you happened to sequence it.

**None of this is somatic calling.** These gates reduce germline contamination in a tumour-only analysis. They do not make your calls somatic. If you have matched normals for every case, call somatic variants properly with a tool built for it and feed the result in as Option A.

### Step 2: decide which genes are candidates

Open `gene_whitelist.txt`. Each line is a gene symbol or a pattern. Lines starting with `#` are comments; blank lines are ignored.

Matching is deliberately not "substring everywhere":

- A line of plain characters (`TP53`, `H3-3A`) is an **exact, case-insensitive symbol match**. This is what stops `RB1` from also pulling in `RB1CC1`.
- A line containing a regex metacharacter (`^H3-`, `^KMT2[ABCD]$`, `^MUC[0-9]`) is treated as a **regex** and matched as a substring, so you can catch a whole family in one line.

**An empty whitelist means "consider every gene"**, which is a reasonable first pass on a new cohort - let recurrence pick the panel, and use `top_n_genes` to keep it readable.

Open `gene_blacklist.txt` the same way. Anything matching a blacklist line is dropped even if it also matched the whitelist - **the blacklist always wins**. The shipped defaults cover the genes that recur for reasons that aren't biology: very long genes (`TTN`, `MUC16`, `OBSCN`) where variant count tracks gene length, polymorphic families (`HLA-*`, olfactory receptors), and non-genes SnpEff happily annotates (`LINC*`, `MIR*`, pseudogenes). They're a starting point. Delete any line whose gene you actually care about.

### Step 3: run Section 1

It reports how many variant rows survived each stage, how many genes were dropped by the whitelist, by the blacklist, and by the recurrence cut, which consequence terms weren't in your class map, and which whitelist patterns **matched nothing at all** - almost always a typo, though on a broad panel most patterns legitimately match nothing in any one cohort.

Then it writes `gene_list.csv`, one row per candidate gene, with `include = TRUE` and a `matched_whitelist_pattern` column showing which line pulled each gene in - useful for spotting an overly broad pattern.

**If `gene_list.csv` already exists, Section 1 will not overwrite it.** It writes `gene_list_draft.csv` instead, so a re-run to pick up whitelist edits never silently destroys manual work. Diff the two and merge whatever you want by hand.

### Step 4: edit the CSV

Open `gene_list.csv` in Excel, a text editor, or R. You can:

- Set `include` to `FALSE` to drop a gene without deleting it - the row stays as a record of "considered, rejected", and the `notes` column is where you say why.
- Delete a row entirely.
- Add a row for any gene you want, whether or not a pattern caught it, **including one the blacklist would have dropped.** Section 2 reads the CSV, not the filters, so the blacklist does not get a second vote. You only need `include` and `gene` - the recurrence columns are looked up and written back into the CSV for you on the next Section 2 run. Only genuinely blank cells are filled; anything you typed, `notes` included, is left alone.
- Reorder rows freely - row order in the CSV doesn't matter, since genes are drawn by recurrence (or by `gene_order` if you set it).

If a gene in your CSV has no plottable variant anywhere in the table - a typo, a symbol your annotation uses differently, or a gene that genuinely isn't hit - it's **reported by name and skipped**, rather than drawn as an empty row.

### Step 5: run Section 2

This reads `gene_list.csv`, takes every row with `include = TRUE`, and re-derives each gene's alterations from the variant table. The CSV stores only metadata, never the variants themselves - that's what makes a hand-edit actually take effect.

For each gene × sample cell it collects every alteration class seen there, orders them most-severe-first by `precedence`, and joins them with `;` - the format `oncoPrint()` expects for co-occurring alterations. Classes with a smaller `height` are drawn on top of taller ones, so a gene that's both frameshifted and missense in the same sample shows both.

Samples with no alteration at all still get a column. If you supply `sample_annotation.csv`, that file sets which samples exist and therefore the denominator, so a sequenced-but-unaltered sample isn't quietly dropped. It does **not** set the left-to-right order: that comes from `oncoPrint()`'s own memo-sort unless you set `sample_order`. Any sample marked as a germline normal by `germline_type_value` is excluded from the figure and from every percentage - it's a reference, not a case.

## Note file format (`gene_whitelist.txt` / `gene_blacklist.txt`)

- One gene symbol or pattern per line.
- A line of plain characters is an exact, case-insensitive symbol match.
- A line containing any of `^ $ . * + ? [ ] ( ) | \ { }` is treated as a regex, matched case-insensitively as a substring.
- Exact symbols and regex patterns can be mixed freely in the same file.
- `#` starts a comment; blank lines are ignored.
- Blacklist beats whitelist, always. The only way to force a blacklisted gene into a figure is to add it by hand as a row in `gene_list.csv`.
- An empty whitelist means "every gene is a candidate".

## `consequence_classes.csv` columns

| Column | Meaning |
|---|---|
| `consequence` | The SnpEff term, exactly as it appears in the `ANN` field (e.g. `stop_gained`). Matched case-insensitively. |
| `class` | The category it collapses into. This is what appears in the legend. Several consequences mapping to one class is the normal case. |
| `plot` | `TRUE`/`FALSE` - whether this class is drawn. A `FALSE` row still documents the mapping; it just doesn't reach the figure. |
| `precedence` | Lower = more severe. Decides which class is listed first in a cell with several alterations. Taken per class as the minimum across its terms. |
| `color` | Hex colour for the class. Taken from the first row of that class. |
| `height` | Fraction of the cell height the bar fills, 0-1. Taller bars are drawn first, so **every plotted class needs a distinct height** - two classes at the same height draw identical rectangles and the second one paints the first out completely, leaving the figure disagreeing with `oncoprint_matrix.csv`. Section 2 warns if it finds a tie. |
| `notes` | Free text, entirely yours. Never read by the script. |

To draw something the defaults don't - structural variants, CNVs from a separate caller - add rows with your own `consequence` values and a new `class`, and make sure your variant table uses those exact terms in its consequence column. Nothing about the script is SnpEff-specific past Section 0.

## `gene_list.csv` columns

| Column | Meaning |
|---|---|
| `include` | `TRUE`/`FALSE` - whether this gene is a row in the figure. Along with `gene`, the only column Section 2 actually requires. |
| `gene` | Gene symbol, matched case-insensitively against your variant table. |
| `median_vaf` | Present only when `cols$ao` and `cols$ro` are set. Median within-sample allele fraction across this gene's variants: near 0.50 or 1.00 suggests inherited. Advisory - nothing filters on it. |
| `n_samples` | How many samples have a plottable alteration in this gene. Auto-filled if blank. |
| `pct_samples` | The same as a percentage of all samples. Auto-filled if blank. |
| `n_variants` | Total plottable variant rows for this gene. Auto-filled if blank. |
| `classes_observed` | Which alteration classes were seen, `;`-separated. Auto-filled if blank. |
| `matched_whitelist_pattern` | Which whitelist line pulled this gene in. Empty for genes you added by hand. |
| `notes` | Free text, entirely yours. Never touched by the script. |

## Input table format

One row per variant, per gene, per consequence. Column *names* are yours - map them in `cols`:

| `cols` entry | Required? | What it holds |
|---|---|---|
| `sample` | yes | Sample identifier. Defines the columns of the figure. |
| `gene` | yes | Gene symbol. Defines the rows. |
| `consequence` | yes | A term matched against `consequence_classes.csv`. |
| `impact` | no | `HIGH`/`MODERATE`/`LOW`/`MODIFIER`. If absent, `impact_whitelist` is skipped entirely and filtering falls to the class map alone. |
| `protein` | no | HGVS protein change, carried through to `output/<name>_long.csv`. |
| `chrom`, `pos` | no | Carried through to the long output. |

`.csv`, `.tsv`, `.csv.gz`, `.rds` and `.qs` are all read without being told which is which.

## Configuration reference

All in the CONFIG chunk at the top of `Build_Oncoprint.rmd`:

| Variable | Meaning |
|---|---|
| `base_dir` | Absolute path to wherever you cloned this repo. **Change this first** - it's the one thing the script can't infer. |
| `out_dir` | Where output lands. Defaults to `base_dir/output`; point it at a project folder if you want the figures to live there. |
| `gene_csv` | Path to the working gene CSV. Defaults to `base_dir/gene_list.csv`. |
| `sample_annotation_csv` | Optional per-sample metadata for the annotation bars. `NA` draws none. |
| `variant_table` | Your long annotated variant table. Defaults to what Section 0 writes. |
| `cols` | Which columns of *your* table mean what. See [Input table format](#input-table-format). |
| `vcf_dir` | Folder of SnpEff VCFs for Section 0. `NA` skips it. |
| `min_qual` | QUAL floor when reading VCFs. Default 30. Skipped, with a message, for a file whose `QUAL` is non-numeric (`.`). |
| `filename_fields` | Named regexes pulling sample metadata out of each VCF's file name. |
| `impact_whitelist` | SnpEff impact levels to keep. `NULL` keeps all. Default `HIGH` + `MODERATE`. |
| `min_samples_mutated` | Drop genes altered in fewer than this many samples. Default 2; 1 keeps everything. |
| `top_n_genes` | Cap on panel size, by recurrence then variant count. `NA` for no cap. Default 40. |
| `germline_type_value` | Value of `cols$type` marking a sample as a normal, e.g. `"Germline"`. `NA` disables gate 1. Normals are excluded from the figure and from every percentage. See [Step 1b](#step-1b-filter-out-inherited-variants). |
| `germline_scope` | `"matched"` subtracts a normal from its own donor only (needs `cols$donor`); `"cohort"` treats all normals as a panel. Default `"matched"`. |
| `germline_min_normals` | Under `"cohort"` scope, how many normals must carry a variant before it is dropped. Default 1. |
| `max_pop_af` | Population-frequency ceiling for gate 2. Default `0.01`. Needs `cols$pop_af`; a column named exactly `AF` is refused. |
| `out_name` | Basename for the written files. |
| `plot_title` | Title drawn above the figure. |
| `sample_order` / `gene_order` | Explicit orderings. `NULL` = memo-sort columns, recurrence-sort rows. |
| `show_pct` | Percent-of-samples column on the left. |
| `background_col` | Colour of a "no alteration" cell. |
| `cell_w` / `cell_h` | Millimetres per sample column / gene row. |
| `fig_w` / `fig_h` | Output size in inches. `NA` derives them from the cell sizes and the matrix shape. |
| `write_png` / `png_dpi` | Whether to also write a PNG, and at what resolution. |
| `annotation_colors` | Named vectors fixing the colours of annotation bars. Anything unnamed gets an automatic palette. |

## Output structure

```
output/
├── variants_annotated.csv.gz   # Section 0 only: the parsed long variant table
├── oncoprint.pdf               # the figure, vector
├── oncoprint.png               # the figure, raster
├── oncoprint_matrix.csv        # gene x sample, cells are ';'-joined classes
└── oncoprint_long.csv          # every variant that made it into the figure
```

`oncoprint_matrix.csv` and `oncoprint_long.csv` are there so the figure is auditable - you can check any cell against the variants behind it without re-running anything, and `_long.csv` is what you'd attach as a supplementary table.

## Troubleshooting

- **`Column(s) not in the variant table`** - the names in `cols` don't match your table. The error prints every column your table actually has; copy from that list.
- **`No variants survived the impact/class filters`** - either your table has no `HIGH`/`MODERATE` rows (check the `impact` column really holds SnpEff impact levels, not consequence terms), or every consequence in it maps to a `plot=FALSE` class. Set `impact_whitelist <- NULL` and re-run Section 1 to see the stage-by-stage counts.
- **A gene you expected is missing** - check the Section 1 report in order: was it dropped by the whitelist (no pattern matched), by the blacklist (an overly broad pattern - `^MUC[0-9]` and friends catch more than you'd think), or by `min_samples_mutated`? If it survived to `gene_list.csv` but isn't in the figure, check `include` is `TRUE` and that Section 2 didn't report it as having no plottable variant.
- **A whitelist pattern "matched nothing"** - reported explicitly by Section 1. Usually a typo, a symbol alias (`H3F3A` vs `H3-3A` - gene symbols get renamed, and your annotation may use either), or the gene genuinely isn't altered in this cohort.
- **Consequences reported as "not in consequence_classes.csv"** - your annotator emits a term the shipped map doesn't cover. Add a row for it. This is expected if you're using VEP terms rather than SnpEff's, or a different caller.
- **`Reading a .qs file needs the qs package`** - `install.packages("qs")`, or re-export your table as `.csv.gz`.
- **Section 0 is slow / runs out of memory** - it's reading every variant in every VCF before filtering. Pre-filter the VCFs (`bcftools view -i 'QUAL>=30'`) or raise `min_qual`, and run it once rather than per draft - that's why it writes a reusable table.
- **The figure is unreadable / enormous** - lower `top_n_genes`, or set `fig_w`/`fig_h` explicitly rather than letting them derive from `cell_w`/`cell_h`.
- **Colours changed between runs** - annotation bar colours are auto-assigned from a palette in sorted level order, so adding a sample with a new `CancerType` can shift them. Pin them in `annotation_colors` for anything going in a paper.
- **Edited a file but nothing changed** - re-run the CONFIG chunk. R keeps old variable values in memory until the chunk that sets them runs again, so editing the file alone doesn't update a value already loaded into your session.

## Citing ComplexHeatmap / SnpEff / vcfR

If you publish a figure produced by this tool, cite the underlying software, not just this script:

**ComplexHeatmap (draws the oncoprint):**
> Gu Z, Eils R, Schlesner M. Complex heatmaps reveal patterns and correlations in multidimensional genomic data. *Bioinformatics*. 2016;32(18):2847-2849. doi:[10.1093/bioinformatics/btw313](https://doi.org/10.1093/bioinformatics/btw313)

> Gu Z. Complex heatmap visualization. *iMeta*. 2022;1(3):e43. doi:[10.1002/imt2.43](https://doi.org/10.1002/imt2.43)

**SnpEff (produces the annotations Section 0 parses):**
> Cingolani P, Platts A, Wang LL, et al. A program for annotating and predicting the effects of single nucleotide polymorphisms, SnpEff. *Fly*. 2012;6(2):80-92. doi:[10.4161/fly.19695](https://doi.org/10.4161/fly.19695)

**vcfR (reads the VCFs):**
> Knaus BJ, Grünwald NJ. vcfR: a package to manipulate and visualize variant call format data in R. *Molecular Ecology Resources*. 2017;17(1):44-53. doi:[10.1111/1755-0998.12549](https://doi.org/10.1111/1755-0998.12549)

For the exact, current citation of any of these as installed on your machine, run `citation("ComplexHeatmap")` or `citation("vcfR")` in R - these change between versions, so don't copy a citation from anywhere (including here) without checking that first.

Also cite your **variant caller** (FreeBayes, GATK, Mutect2 - whatever produced the VCFs) and your **reference genome build**. Neither is visible to this tool, and both matter more to the figure than anything in this repo does.

## Citing this repository

This repo ships a [`CITATION.cff`](CITATION.cff) file, which is GitHub's standard machine-readable citation format - GitHub shows a **"Cite this repository"** button on the repo's main page (top right, near "About"), which generates APA and BibTeX for you and stays in sync if you edit `CITATION.cff`. That's the recommended way to get a citation - it's always current with whatever's in the file.

For quick reference, the same details in a few common formats:

**Plain text:**
> Balagopal, N. K. (2026). *Oncoprint Builder* [Computer software]. GitHub. https://github.com/navinvent/Oncoprint_Builder

**APA:**
> Balagopal, N. K. (2026). *Oncoprint Builder* (Version 1.1.0) [Computer software]. GitHub. https://github.com/navinvent/Oncoprint_Builder. ORCID: [0000-0002-3012-2469](https://orcid.org/0000-0002-3012-2469)

**BibTeX** (`@software` - if your BibTeX style doesn't support the `software` entry type, e.g. plain `bibtex` rather than `biblatex`, use `@misc` instead, shown second):

```bibtex
@software{balagopal_2026_oncoprint_builder,
  author  = {Balagopal, Navin Kumar},
  title   = {{Oncoprint Builder}},
  year    = {2026},
  url     = {https://github.com/navinvent/Oncoprint_Builder},
  version = {1.0.0}
}
```

```bibtex
@misc{balagopal_2026_oncoprint_builder,
  author       = {Balagopal, Navin Kumar},
  title        = {{Oncoprint Builder}},
  year         = {2026},
  howpublished = {\url{https://github.com/navinvent/Oncoprint_Builder}},
  note         = {Version 1.0.0}
}
```

**One more option worth knowing about:** if you want this tool to have a permanent, versioned citation with its own DOI (common for research software people actually depend on, and required by some journals for "software availability" statements), you can archive a GitHub release through [Zenodo](https://zenodo.org) - link your GitHub account, flip the switch for this repo, then cut a release (a git tag) on GitHub. Zenodo mints a DOI for that exact release automatically and it's citable forever, independent of whether the GitHub repo changes or disappears later. Not set up yet here since there hasn't been a tagged release.

## License

[MIT](LICENSE) - permissive, allows reuse/modification/redistribution including commercially, provided the copyright notice is kept.

## Provenance

Written to fill a gap in an internal radiation-induced glioma whole-exome pipeline, whose VCF import and `ANN=` parsing stages produced a long annotated variant table and then stopped - the oncoprint stage existed only as a placeholder README saying no code existed yet. The SnpEff parsing in Section 0 follows that pipeline's approach (QUAL ≥ 30, chunked `ANN=` expansion, `&`-split consequences, sample metadata derived from file names); everything downstream of the variant table is new. The whitelist/blacklist split, the reviewable intermediate CSV, and the "the CSV is the source of truth, not the filters" rule are carried over from [GMT Collection Builder](https://github.com/navinvent/GMT_Collection_Builder), which solves the same problem for pathway panels.

No internal data, VCFs, or cohort metadata are included in this repo. `gene_list.example.csv` and `sample_annotation.example.csv` are illustrative placeholders showing file shape - the numbers in them are not results.
