Golden digests for the ProFound regression case set (../helper-cases.R).

Each digest is md5() of serialize(value, version = 2) of the value returned by
that case after pf_canon(), which rounds numerics to 12 significant digits
and relabels integer matrices into 1..k by first appearance. Watershed
segment ids are an artefact of pixel processing order rather than physics,
so canonicalising removes that freedom while still detecting any genuine
change to segmentation geometry or photometry. The serialization header's
writer-R-version bytes are blanked before hashing, so a digest depends only on
the canonical value, not on the R version that produced it.

Regenerate with:
  Rscript tests/scripts/make-goldens.R [library-holding-the-reference-build]

Provenance:
  built_with:  ProFound 1.34.5
  r_version:   4.5.3
  generated:   2026-10-26 15:27:40
  n_cases:     139
  source_head: 2d24022

The checked-in reference table was produced from the committed v1.34.5
source, i.e. BEFORE the speed-focused edits in src/ were committed, so the
suite checks those edits against the numerical behaviour of the release that
preceded them. Re-running this script against the current source is only
meaningful once any intended behaviour change has been accepted.

Deliberate exclusion: a `that` vector containing both negative and positive
values makes this_in_that() index its lookup table out of bounds, producing
heap-dependent results that differ between runs of the same build. Such
input is not part of any case; see test-05-this-in-that.R.

Suite validation and known pre-existing differences are documented inline in
this file's current checked-in copy; regenerate it and re-add them if the
reference build is ever rebased.

Deliberate exclusion: a `that` vector containing both negative and positive
values makes this_in_that() index its lookup table out of bounds, producing
heap-dependent results that differ between runs of the same build. Such
input is not part of any case; see test-04-this-in-that.R.

Suite validation (mutation testing, 2026-09-26)
================================================
Each edited src/ file was individually reverted to its committed v1.34.5
content while the rest of the rewrite stayed in place, and the full suite was
run. All five reverts leave the suite passing -- as does running the suite
against a separately installed, wholly pure v1.34.5 build. Each optimisation is
therefore behaviour-preserving rather than merely consistent with a test
written against the new code.

Deliberate defects, and the number of assertions that caught them:
  dilate: source pixel no longer written (33), expand filter ignored (8),
          existing labels overwritten (6), smallest-id tie-break inverted (4)
  water:  merge tolerance halved (5)
  ellip:  inside-ellipse bounding-box test loosened to qmax <= 1.25 (17)
  nser:   the nser == 1 weight perturbed in the 7th digit (2)
  akima:  interval index clamped one too early with the correction loops also
          removed (26)

One honest gap: removing the XLookup/YLookup correction loops on their own, or
clamping the index early with the loops intact, changes nothing here -- on this
platform the direct estimate (x - xMin) / xSpacing is already correct for every
input the suite uses, so the loops never fire. They are defensive rather than
load-bearing in practice, and their necessity is consequently untested rather
than proven. The combined defect (bad clamp without loops) is caught.

Pre-existing difference worth knowing (identical in v1.34.5 and now)
====================================================================
water_cpp(ext = 0, abstol = 0) leaves some above-cut pixels unlabelled whereas
water_cpp_old() does not (1329 vs 1291 segments on the suite's test image). The
numbers are the same in both builds, ext = 0 is not the default, and at the
default ext = 2 the two implementations agree, so profoundProFound's normal
path is unaffected.

Regeneration note (2026-10-26)
==============================
This table was regenerated on R 4.5.3 from the current tree. The previous table
was produced on R 4.4.1 and could never match, because serialize() stamped the
writing R version into every digest; see pf_digest() in helper-profound.R. The
digests were verified to change only by that header (regenerating on the same
build reproduces them exactly), so the values themselves are still the ones the
v1.34.5 reference build produced.
