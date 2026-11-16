# Presentation Cheat Sheet

Keep this open during questions. Everything here is in the notebook or the README.

## The project in one sentence

We compared the COVID virus's spike protein across the Alpha, Delta and Omicron
variants, tested whether the differences sit in one functional region, and
trained models to predict a sample's geographic region from its sequence.

## Numbers worth knowing by heart

| | |
|---|---|
| Samples downloaded | 14,965 |
| Distinct sequences analysed | 1,406 |
| Positions in spike | 1,273 |
| Positions that vary at all | 583 |
| Positions significantly different between variants | 83 |
| Median mutations in the binding region | Alpha 1, Delta 2, Omicron 15 |
| Accuracy predicting variant | 1.000 (circular, a control only) |
| Accuracy predicting region, B.1.1.7 / B.1.617.2 / BA.2 | 0.659 / 0.832 / 0.944 |
| Guessing baseline | 0.500 |
| Highest score any sequence-only model could reach | 0.871 |
| Random positions needed to identify a variant | about 50 |

## The questions you will probably get

**Why isn't this the question on your approval form?**
The original question was answerable without data. Variants are defined by their
spike mutations, so counting them in a published catalogue answers "does Omicron
have more" and a classifier predicting variant is recovering the labelling rule.
Slide 2 has the catalogue counts. We kept the comparison as descriptive work and
built the analysis around questions that were genuinely open.

**Your accuracy is 1.000. Isn't that overfitting?**
No, it is circularity, which is different. The model is not memorising noise, it
is reading the exact mutations used to assign the labels. A majority-class
baseline scores 0.333 and clustering with no labels at all recovers the variants
at an adjusted Rand index of 0.998. We report it as a control, not a result.

**How do you know there is no data leakage?**
Three ways. Identical sequences were removed before splitting, so no sequence
appears in training and test. Feature selection runs inside each
cross-validation fold through a pipeline, never on the full dataset. The test
set is scored once.

**Why balanced accuracy rather than accuracy?**
The classes are uneven. Plain accuracy would reward a model for always
predicting the larger region. Balanced accuracy averages the score across both
classes, so guessing scores 0.500 by construction.

**How did you handle testing 1,273 positions at once?**
Benjamini-Hochberg correction, and we report q-values rather than p-values. We
also report the effect size beside every q-value, because with hundreds of
sequences a difference of a few percent can be significant while meaning little.

**Why only Europe and North America?**
After removing duplicate sequences, too few distinct sequences remained from
Africa and Asia to cross-validate. That reflects where sequencing was done, not
where the virus was. It is a limitation, and it is in the notebook.

**What is the 0.871 ceiling?**
Some identical spike sequences were collected on different continents. No model
that sees only the sequence can assign all of those correctly. We measured that
bound and report it next to the model scores, so the 0.659 is read against 0.871
rather than against a perfect score.

**What does the receptor binding domain actually do?**
It is the part of spike that attaches to human cells. Positions 319 to 541.
Changes there can alter how tightly the virus binds, and can help it avoid
antibodies produced against earlier variants.

**What is the furin cleavage site finding?**
Positions 681 to 685, where the protein must be cut for the virus to fuse with a
cell. Alpha changes it one way and Delta another. Two variants altering the same
five residues by different means points to the same pressure acting twice, not
to shared ancestry.

**How do you know your pipeline is correct?**
We never fed the published mutation catalogue into the analysis, and the
medians it produced, 1, 2 and 15, match published figures. There is also a
validation suite checking the reference length, twelve known residues, the
domain boundaries and agreement with the catalogue for every group used.

**Why did you remove 13,559 of 14,965 samples?**
Mostly duplication. Many samples have identical spike sequences, and keeping
them would place copies of the same sequence in both training and test data. The
full filter chain with before and after counts is in Stage 3.

**What would you do differently?**
Align the DNA rather than the protein, which would place deletions correctly in
one region where protein alignment cannot. Build a phylogenetic tree to test
whether the geographic signal survives accounting for how samples are related.

## If you are asked something you do not know

Say so, then say where the answer is. "That is in Stage 5 of the notebook, I can
open it." Guessing in front of a marker costs more than admitting a gap.

## Do not say

- "The model predicts the variant with 100% accuracy" as if it were an
  achievement. Always add that the task is circular.
- "It proves geography causes the mutations." It does not. Related samples
  cluster together, and we cannot separate that from geography.
- Any number not on this sheet. Check the notebook instead.
