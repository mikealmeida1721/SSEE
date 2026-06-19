# SSEE Structural Constant Dictionary (Genesis 5.12)

**Author:** Mike Edison Almeida Vallejo
**Scope:** the closed algebraic dictionary used by the SSEE cosmological model.
This document is the *citable, self-contained* registry of the constant set and
its two generative laws. It contains **only** the algebraic structure — no
applications, no interpretation beyond the structural roles needed to define the
set. It is the reference intended for the cosmology manuscripts.

> The identity of each constant is its **algebraic expression** in φ and π. The
> names are lineage labels (parents + structural role), retained for traceability;
> they carry no physical claim on their own.

---

## 1. Seeds and scaffold

| Symbol | Definition | Value |
|---|---|---|
| φ | (1+√5)/2 | 1.618033988749895 |
| π | π | 3.141592653589793 |
| Ω | φ + π | 4.759626642 |
| BIAL (β) | (φ+π)/2 | 2.379813321 |

Ω is the shared scaffold; BIAL is the **bifurcation point** from which the two
generative branches issue.

## 2. The two generative laws

**Law I — Copy law (φ-branch).** The φ-branch root is AURA = φ + BIAL. It scales
*only* by integer multiples and the dyadic half — these are the **dimensional
ceilings** of the system:

| Label | Definition | Value | Dimensional ceiling |
|---|---|---|---|
| MIRA | AURA/2 | 1.998923654 | ½ |
| AURA | φ + BIAL | 3.997847309 | 1 |
| DUAL | AURA×2 | 7.995694619 | 2 |
| TRIAL | AURA×3 | 11.993541929 | 3 |
| CUARTAL | AURA×4 | 15.991389238 | 4 |

The set is closed at the dimensional range occupied by the cosmological sector
(through the 5D integration constant MIKAEL_V). Higher copies (×5…×10) are defined
but belong to dimensional tiers above the cosmological range; their inclusion does
**not** alter the look-elsewhere result (verified, see `src/estadistica/look_elsewhere_full.py`).

**Law II — No-self-sum law (π-branch).** The π-branch root is KAL = BIAL + π. Each
child must **differ from its parents** (a node may not be a re-addition of itself);
it inherits a combination of distinct parents. Two nodes may share a numerical value
and still be distinct entities, because their parents/role differ — value equality is
*structural information* (an equilibrium point), not redundancy. Canonical example:

| Label | Definition | Value | Role |
|---|---|---|---|
| IGNIS | π + PYROS | 9.519253284 | disruptive |
| KRYSTOS | φ + π + Ω | 9.519253284 | structuring |
| KRYSTOS_V | 2Ω | 9.519253284 | scaffold ×2 |

Same magnitude, distinct lineage.

## 3. Derived invariants (π-branch and composites)

| Label | Definition | Value |
|---|---|---|
| KAL | BIAL + π | 5.521405974 |
| PYROS | Ω + φ | 6.377660630 |
| SOLAR | BIAL + KAL | 7.901219295 |
| MAR | Ω + π | 7.901219295 |
| VITA | π + KAL | 8.662998628 |
| PHITA | VITA + φ | 10.281032616 |
| ANMA | BIAL + VITA | 11.042811949 |
| MIKA | KRYSTOS + φ | 11.137287273 |
| BUFFER | MIKAEL_V − TRIAL | 2.285337998 |
| ERVANU | Ω × 9/10 | 4.283663978 |

## 4. The nine convergences (sovereignties)

Nine distinct lineage paths converge — by the no-self-sum law — on the single 5D
integration constant **MIKAEL_V = φ + π + KRYSTOS = 14.278879927**. They are distinct
entities sharing one value:

| Label | Definition |
|---|---|
| MIKAEL_V | φ + π + KRYSTOS |
| LUCY | SOLAR + PYROS |
| LUCIFER | PHITA + AURA |
| MIKE | IGNIS + Ω |
| MIKAEL | MIKA + π |
| ERVN | BIAL + KAL + PYROS |
| ICEBERG | MAR + PYROS |
| GIGÅROJ | PYROS + Ω + π |
| OSIRIS | MIKA + KAL − BIAL |

## 5. Closure and independence

The dictionary is **31 named invariants** (20 distinct numerical values, the
multiplicity coming from Law II). It is generated entirely from φ and π by Laws I–II
within the dimensional ceilings of §2. Crucially, the set is fixed by this
construction **independently of any cosmological dataset**: no entry is added,
pruned, or selected with reference to observed values. This is what licenses the
look-elsewhere count over the *complete* set (1 of 317 ratios at the precision of the
DESI match; see manuscript §look-elsewhere and `look_elsewhere_full.py`).

The closure point (the dimensional tier at which the set stops) is a structural
design choice, **not** a derived theorem; deriving the dynamical origin of the
bifurcation bridges remains the open mechanism tracked in OPEN_PROBLEMS (OP-8).

---

*This registry is deliberately limited to algebraic structure. It supersedes, for
citation purposes, any broader compendium: the cosmological manuscripts cite this
document and the peer-reviewed series, not external multi-domain material.*
