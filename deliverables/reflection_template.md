# Individual Reflection — Computing Bioinformatics Group Project

**Name:**  ____________________    **Admission number:**  ____________________

**Project:** SARS-CoV-2 spike protein variation across Alpha, Delta and Omicron

> **This document is a scaffold, not a submission.** The reflection is assessed
> individually and asks what *you* did. The prompts and the factual summary
> below are here so you are not starting from a blank page. Replace every
> prompt with your own account before submitting. Keep it to one page.

---

## 1. My role in the group project

*What were you responsible for? Be specific about tasks, not job titles. Which
stage or stages did you own? What did you produce that the group used? If you
coordinated rather than coded, say what coordinating actually involved.*

____________________________________________________________________

## 2. Challenges experienced

*What actually went wrong, and what did you do about it? Specific beats
general. Below are real difficulties this project ran into — use any that you
personally dealt with, and add your own.*

- The data source did not behave as documented. The NCBI download endpoint
  returned metadata but no sequences for specific accessions, and the
  whole-lineage alternative was about 641 MB, so the sequences had to be
  fetched a different way.
- Most Alpha genomes are flagged `PARTIAL` even when the spike gene is
  complete, so the obvious quality filter would have thrown away most of one
  variant.
- Around 7,350 of the records turned out to be duplicates of other records at
  the spike level, which would have leaked between training and test data.
- The first attempt at deduplication counted sequences as distinct when they
  differed only in where their unresolved positions fell, which inflated the
  apparent diversity.
- Geographic region turned out to be confounded first with variant, and then,
  less obviously, with Pango sublineage inside each variant.

____________________________________________________________________

## 3. Achievements gained

*What can you now do that you could not before? Name specific skills, and
point to something concrete in the notebook that shows it.*

Candidate areas, if they match what you did:

- Working with a real biological database API, including handling the case
  where its documented behaviour and actual behaviour differ.
- Sequence alignment and coordinate systems: why a mutation is called S:501
  and what has to be true for that number to mean anything.
- Quality control reasoning: choosing filters based on what they do to the
  biology rather than applying defaults.
- Statistical testing with multiple-comparison correction, and reporting
  effect size separately from significance.
- Designing a machine learning experiment that avoids leakage, and recognising
  when a high accuracy score is measuring the wrong thing.
- Reporting an honest negative or modest result instead of overstating it.

____________________________________________________________________

## 4. What I would do differently

*One honest paragraph. This usually scores better than claiming everything
went well.*

____________________________________________________________________

---

### Factual summary of the project, for reference

The project compared spike protein variation between Alpha, Delta and Omicron,
asked whether the variant-distinguishing signal is confined to the
receptor-binding domain, and tested whether an isolate's geographic region can
be predicted from its spike sequence.

14,965 isolate records were retrieved from NCBI and filtered to 1,406 distinct
spike mutation profiles. Mutation frequency differs significantly between
variants at 83 positions after Benjamini-Hochberg correction. Median RBD
mutations per isolate are 1 in Alpha, 2 in Delta and 15 in Omicron. Masking
the RBD does not reduce variant classification accuracy at all, and about
fifty arbitrary varying positions are enough to identify a variant, so the
signal is redundant rather than localised. Geographic region is predictable
above baseline within a single Pango lineage, replicated in B.1.1.7,
B.1.617.2 and BA.2.
