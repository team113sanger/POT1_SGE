# POT1 molecular dynamics simulations

All-atom molecular dynamics simulations of wild-type POT1 and six clinically observed
missense variants bound to telomeric single-stranded DNA — one microsecond per system,
seven systems.

These are exploratory simulations, shared alongside the paper. They are not part of its
analysis, and no result in the paper depends on them. They are provided because the
question they ask is a reasonable one and the groundwork may be useful to someone pursuing
it properly.

Companion to Obolenski *et al.*, *A Functional Spectrum of POT1 Alleles Dictates
Genotype-Phenotype Associations*.

---

## The question

Saturation genome editing (SGE) places every POT1 missense variant on a quantitative
functional scale and separates the depleted variants into strong and weak classes. That
classification is empirical: it says how a variant behaves in the assay, not why.

These simulations asked whether the distinction has a visible structural correlate — whether
variants in different functional classes perturb the dynamics of the POT1–ssDNA complex in
different ways. Six variants were chosen to span the classes and to span two structural
situations: three that contact the telomeric DNA directly, and two that sit away from it at
positions important for the fold.

---

## The simulations

![Normal mode and RMSF analysis for six POT1 missense variants, with insets showing hydrogen bonding and energy decomposition values relative to wild type](figures/protein_dynamics.svg)

> Individual structures show the first normal mode (NM1), with magnitude of motion for NM1
> depicted by the width of the tube, for each key missense variant. Each structure is
> coloured as a heatmap based on the change in RMSF (ΔRMSF) with respect to the wild-type
> RMSF. Insets highlight significant changes in hydrogen bonding interactions (dashed lines)
> with the telomeric DNA with respect to wild-type, and changes in the intermolecular
> interactions as calculated from energy decomposition analysis (EDA) with respect to the
> corresponding mutation of each variant (individual values under each residue; positive
> values indicate interactions for wild-type are favoured and negative values indicate the
> variant is favoured).

The wild-type system is the reference throughout: it sets the baseline RMSF from which
ΔRMSF is measured and the hydrogen bonding and EDA values are compared.

### Systems

Built from PDB [1XJV](https://www.rcsb.org/structure/1XJV) — POT1 OB1 and OB2 domains bound
to telomeric ssDNA.

| System | SGE functional classification | ClinVar (at extraction) | Position relative to ssDNA |
|---|---|---|---|
| Wild type | — | — | — |
| p.Y36C | Weak depleted | Uncertain significance | Distant from DNA |
| p.D51N | Strong depleted | Uncertain significance | Distant from DNA; site important for protein fold |
| p.Y89C | Weak depleted | Conflicting interpretation | Direct DNA interaction |
| p.Q94E | Weak depleted | Uncertain significance | Direct DNA interaction |
| p.Q94R | Weak depleted | Uncertain significance | Direct DNA interaction |
| p.R137H | Strong depleted | Uncertain significance | Distant from DNA; site important for protein fold |

The SGE classifications come from the editing screen and stand on their own. They are in
`supplementary_data/` in this repository.

### Methods as performed

| | |
|---|---|
| Engine | AMBER20, `pmemd.cuda` |
| Force field | ff19SB |
| Starting model | PDB 1XJV; Coot and MolProbity used to alleviate bad clashes |
| Missing residues | Modeller interface in UCSF Chimera |
| Variant construction | Dunbrack rotamer library in Chimera |
| Solvation | TIP3P water, cubic box extending 12 Å from complex surface |
| Counterions | Potassium, to neutralise system charge |
| Ensemble | NPT, Langevin thermostat |
| Non-bonded cutoff | 9 Å, smooth PME |
| Production length | 1 μs per system |
| Independent replicas | 1 per system |
| Trajectory analysis | `cpptraj` (RMSD, RMSF, hydrogen bonds) |
| Normal mode analysis | `cpptraj` and ProDy (VMD) |
| Energy decomposition | AMBER-EDA |
| Binding enthalpies | MM-GBSA via `MMPBSA.py` (AmberTools) |
| Downstream analysis | R (`data.table`, `tidyverse`, `abind`) |
| Visualisation | VMD, UCSF Chimera, ChimeraX, gnuplot, matplotlib, ggplot2 |

---

## What can be read from this

The useful part of this dataset is the setup, not the comparisons. Read in the right
register it is informative; read as measurement it is not.

**Each system was simulated once.** With one trajectory per system there is no estimate of
within-system variation, so a difference between two systems cannot be separated from the
variation either would show on being run again. Everything comparative here — ΔRMSF, normal
mode amplitude, hydrogen bond occupancy, EDA energies — is therefore an observation to be
tested, not a quantity that has been measured. This is the constraint that governs all the
others.

**The numbers in the insets have no dispersion.** They are single-trajectory values, so no
uncertainty can be attached to them. EDA-derived interaction energies are in any case more
reliable as relative than absolute quantities, and are strongest alongside experimental
validation.

**Convergence was not characterised.** No RMSD-against-time analysis was performed, so
whether 1 μs is long enough for these systems to equilibrate is not established.

**Starting geometries were modelled.** Mutant side chains were built with the Dunbrack
rotamer library and clashes relieved. Rotamer choice matters at these positions, especially
the three that contact DNA, so some of what is seen may follow from the initial model rather
than from the substitution.

**The analysis sits at the level of global motion.** Normal modes and whole-chain RMSF
describe how the complex moves overall. For substitutions like these the more informative
level is local — what happens to a specific contact when a residue changes. The clearest
example is p.Y89C, where the π–π stacking between Y89 and the telomeric guanine is lost
outright. That is a concrete, tractable question this dataset is positioned to address but
does not yet answer.

**What the figure does show,** taken as orientation rather than as result, is where the six
substitutions sit relative to the DNA, which contacts each one is near, and which regions of
the complex move most in the first normal mode. That is a reasonable starting point for
designing a properly powered study.

---

## Further work

The obvious next step is replication. Most of the cost of this kind of study is in building,
parameterising and equilibrating the systems, and that work is done.

- **Multiple independent replicas per system**, initialised from different random velocity
  seeds — commonly more than five — so that between-system differences can be separated from
  within-system variation. This single change converts every observation above into a
  testable quantity.
- **Convergence assessed and shown**, with RMSD-against-time and related diagnostics, rather
  than inferred from trajectory length.
- **Local interaction analysis** at the substituted positions, the Y89–guanine stacking being
  the first case to take: occupancy and geometry of that interaction across replicas, against
  wild type.
- **Sensitivity to the starting model**, by sampling rotamer choice at the substituted
  positions rather than committing to a single built geometry.
- **Energetics with uncertainty**, reporting hydrogen bonding and EDA-derived interaction
  energies across replicas and framing them as relative quantities.

Whether the strong and weak depleted classes have a structural basis is, on this evidence,
an open question. It is a good question, and this is a reasonable place to start on it.

---

## Citation

The molecular dynamics component does not appear in the paper. If referring to this work,
please cite the primary paper:

> Obolenski S., Waters A.J., *et al.* A Functional Spectrum of POT1 Alleles Dictates
> Genotype-Phenotype Associations.
> <!-- EDIT: add journal, year, DOI on acceptance -->

---

## Contact

<!-- EDIT: confirm contact for queries about the simulations -->
Please open an issue on this repository.
