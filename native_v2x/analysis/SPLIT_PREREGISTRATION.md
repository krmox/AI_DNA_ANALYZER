# V2.x development / holdout protocol (written BEFORE any feature table or label join was inspected)

* Region: HG002 chr20, 300x, GRCh38, NISTv4.2.1 truth, chr20_highconf.bed, hap.py 0.3.15 (unchanged).
* Split unit: 2,000,000-bp genomic blocks, block = floor(pos0 / 2,000,000) (32 blocks, 0..32).
  * DEVELOPMENT = even blocks (0,2,4,...).   HOLDOUT = odd blocks (1,3,5,...).
* Everything that is fitted or chosen (feature subset, classifier coefficients, uncertainty band, adaptive
  MAX_READS rule, context window, thresholds) is fitted / chosen ONLY on DEVELOPMENT loci.
  Cross-validation for model selection is by block (leave-blocks-out) inside DEVELOPMENT.
* HOLDOUT is scored once per frozen configuration ("finalist"); a finalist is frozen (coefficients written to
  a model file, its SHA-256 recorded) before its holdout number is computed.  No config is edited after seeing
  its holdout number; failed finalists are reported, not replaced silently.
* Blocks are contiguous 2 Mb regions, larger than any read / context window (<= 25 bp) used, so no read or
  context window straddles the two sets except at block borders (border effect: < 100 bp per border).
* Uncertainty: block bootstrap (resample whole 2 Mb blocks) for all deltas vs V2 control.
* Full-chr20 numbers are reported as ENGINEERING benchmarks (dev+holdout mixed) and are labelled as such;
  only the HOLDOUT rows are independent validation.  Frozen constants (Binomial 7.0, router cutoff
  5.411872376933351, PB 10.5, Method C eps 0.01) and the V2 control (MAX_READS=48) are never modified.
* Limitation stated in advance: dev and holdout come from ONE individual and ONE chromosome; that is
  within-sample generalisation across regions, not across samples or chromosomes.
