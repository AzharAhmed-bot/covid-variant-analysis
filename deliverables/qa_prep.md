# Question Preparation

Short answers you can actually say out loud. Each one is two or three sentences.
If you need more, the detail is in the notebook stage named at the end.

---

## Near certain

**Why isn't this the question on your approval form?**
Because the original question could be answered from a lookup table. Variants
are defined by their spike mutations, so "does Omicron have more" is reading the
catalogue, and a classifier predicting variant is learning the labelling rule.
We kept both as checks and built the project on questions that were actually
open. (Slide 2)

**Your accuracy is 1.000. Isn't that overfitting?**
No, it's circularity, which is different. Overfitting is memorising noise. Here
the model is reading the exact mutations used to assign the labels. A baseline
gets 0.333 and clustering with no labels at all recovers the variants at 0.998.
(Stage 7)

**So what did you actually find?**
Three things. The differences sit in functional regions of the protein. The
variant signal is repeated across the whole protein rather than held in the
binding domain. And sequence carries partial information about geography, which
survives once you control for lineage. (Slide 12)

**Is 1,406 sequences enough?**
For the variant comparison, comfortably. For the geography task it's thin: the
three lineage tests use between 229 and 529 sequences each. That's why we report
cross-validation spread and a held-out test score rather than a single number.
(Stage 7)

**Why only Europe and North America?**
After removing duplicates, Africa and Asia had too few distinct sequences to
cross-validate. That reflects where sequencing capacity was, not where the virus
was. It's a limitation, and it's on the limits slide. (Slide 11)

---

## Methods

**Why balanced accuracy instead of accuracy?**
The classes are uneven. Plain accuracy would reward a model for always
predicting the larger region. Balanced accuracy averages recall across classes,
so guessing sits at 0.500 by construction.

**How did you handle multiple testing?**
Benjamini-Hochberg across all tested positions, and we report q-values rather
than p-values. We also report the effect size beside every q-value, because with
over a thousand sequences a difference of a few percent can be significant while
meaning very little. (Stage 5)

**Your p-value is ten to the minus 258. Isn't that absurd?**
It's what you get from a large sample and near-complete separation, and on its
own it isn't informative. That's why we report the rank-biserial correlation
too, which is minus one: every Omicron isolate carries more RBD mutations than
every Alpha one. The effect size is the meaningful number.

**How is the data represented for the model?**
A 1,406 by 1,273 binary matrix. One row per sequence, one column per position,
one if the residue differs from the reference, zero if it matches. Unresolved
positions are encoded zero so missing data can't look like signal.

**How do you know there's no leakage?**
Three ways. Identical mutation profiles are removed before splitting, so no
sequence appears in both halves. Feature selection runs inside each
cross-validation fold through a pipeline, never on the full data. The test set
is scored once. (Stage 7)

**Why both logistic regression and random forest?**
As a check on each other. If a linear model and a tree ensemble disagree
sharply, something is usually wrong with the features. Here they agree closely,
which is reassuring.

**Did you tune hyperparameters?**
No. We fixed them in advance and left them alone. That means the numbers aren't
optimised, but it removes another route for test information to leak into
training, which mattered more to us.

**Why remove duplicates? Isn't that throwing away real data?**
At the sequence level they carry no extra information, and keeping them would
put copies of one sequence in both training and test. We do report what they
tell us: 77 duplicate groups span more than one region, which is where the
0.871 ceiling comes from. (Stage 3)

---

## Harder, probing questions

**What is the 0.871 ceiling and how did you get it?**
Some identical spike sequences were collected on different continents. For each
distinct sequence the best any model can do is pick its commonest region;
everything else is necessarily wrong. Summing that over the dataset gives 0.871.
It's a property of the data, not of our model.

**Couldn't the geographic signal just be collection date?**
We checked. Within each lineage the two regions were sampled over essentially
the same period, with median dates within a few weeks. For Alpha both peak in
the second quarter of 2021. It's on the slide. (Stage 7)

**Have you really controlled for phylogeny?**
Only partly, and we say so. Holding the Pango lineage fixed removes the
sublineage shortcut, but a named lineage still contains finer structure.
Sequences from one region may be related for transmission reasons rather than
geographic ones. Settling it needs a whole-genome tree with ancestral state
reconstruction. (Slide 11)

**Why is the signal peptide not on your domain chart?**
It actually scores highest per residue, about 69 per 100. We left it off because
it's a 13-residue stretch cleaved from the mature protein, sitting at the start
of the gene where deposited sequences are often truncated. In our pipeline a
truncated start reads as deletions, so we treat it as a technical artefact. It's
in the notebook with that caveat.

**How do you know your alignment is correct?**
We validated the reference against twelve known residues and against every
mutation the CoVariants catalogue annotates for the clades we use, with zero
mismatches. We also know where it fails: it misplaces the BA.1 deletion around
142 to 145, because BLOSUM scores the substituted residue better at position 142
than at 145. We declare those windows. (Stage 3)

---

## Where the honest answer is "we didn't"

Don't bluff these. Saying so is fine.

- **Did you build a phylogenetic tree?** No. That's the main thing we'd add next.
- **Did you use the nucleotide sequences?** No, protein only. Aligning codons
  would fix the deletion placement problem.
- **Did you model insertions?** No. They're counted but excluded from features,
  to keep position numbering matching the literature.
- **Did you test other classifiers, or deep learning?** No. With 1,273 binary
  features and about 1,400 rows, more capacity wasn't the bottleneck.

---

## If you're asked how you built it

Answer according to your department's policy on tool use, honestly. The thing
that actually protects you in a viva is being able to explain any decision in
your own words: why completeness wasn't used as a filter, why deduplication keys
on the mutation profile, why the first geography result was wrong. If you can do
that, you understand the project. Work through this sheet until you can.

---

## If you don't know

Say so, then say where it lives. "That's in Stage 5 of the notebook, I can open
it." A marker will take that over a guess every time.
