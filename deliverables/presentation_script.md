# Presentation Script

Single presenter. Roughly 11 to 12 minutes at a measured pace.

Say it in your own words. Keep the order and keep the numbers exact.

---

## Slide 1: Title

Good morning. I'm Azhar, presenting our group project on the SARS-CoV-2 spike
protein.

Those red parts you can see are the spike protein. It's what the virus uses to
get into human cells, and it's what changes most between variants.

Our finding is in the title. If you give us a virus sequence, telling you which
variant it is turns out to be easy. Telling you where in the world it came from
is much harder. I'll explain both halves of that.

## Slide 2: Why the question changed

First, something I have to be upfront about. The question we got approved is not
the question we answered, and here's why.

Variants like Alpha and Delta aren't discovered by looking at data. They are
defined by a list of mutations. A sequence is called Omicron because it carries
the Omicron mutations.

So when our original question asked whether Omicron has more mutations in the
binding domain, the answer is already in the catalogue. Alpha has one, Delta has
two, Omicron has fifteen to seventeen. We didn't need to collect anything to
know that.

The second half had the same problem. Training a model to predict the variant
means the model just learns the rule that was used to apply the labels in the
first place.

So we kept that comparison, but only as a check that our pipeline works. The
real project is built on questions where we genuinely did not know the answer.

## Slide 3: The three questions

We settled on three.

One: where along the protein do the variants actually differ? Not just at the
catalogue positions, but across all 1,273 of them.

Two: is that difference concentrated in the binding domain, which is the part
everyone focuses on? We tested that by deleting the binding domain from the data
and seeing how much accuracy we lost.

Three, and this is the one that matters: can a model tell which part of the world
a sample came from, using only the sequence?

That third one is the honest question, because where a sample was collected has
nothing to do with how a lineage is defined. So the model has to find real
signal. It could have found nothing, and we agreed beforehand that we would
report that if it happened.

## Slide 4: The data

On the data. We pulled 14,965 records from NCBI.

After requiring a sequence and a usable collection date, and dropping sequences
with too many unresolved positions, we were left with 8,756.

Then the big cut, down to 1,406. That's because most samples of the same variant
are identical in spike. Thousands of records collapse to a few hundred genuinely
different sequences.

We had to remove those duplicates, because if the same sequence appears in both
the training and the test set, a model can score well just by remembering it.

And that was the real limit on this project. Not how much we could download, but
how many distinct sequences exist. We went from 2,600 records to 15,000 and
roughly doubled the distinct sequences, not six times.

## Slide 5: Mutation counts

This is the check I mentioned.

On the left, total mutations per sample. On the right, only those in the binding
domain. Each dot is one virus sample.

The medians are one for Alpha, two for Delta, and fifteen for Omicron.

Now, those are the published numbers. But we never gave the catalogue to our
code. It worked them out from the sequences alone and landed on the same values.

That matters, because if this slide had come out wrong, nothing after it would
be worth showing you. The statistics are at the bottom. The difference between
the three groups is significant, and Omicron separates completely from the other
two.

## Slide 6: Where the changes fall

Now where along the protein those changes sit. The axis runs across all 1,273
positions, one row per variant.

Two things to notice. The changes aren't spread out evenly. Each variant has a
small set of positions that nearly every sample has changed, and long stretches
that never change at all.

And Omicron has that dense band through the shaded region, which is the binding
domain. Alpha and Delta don't.

We tested every position separately, and corrected for the fact that we ran over
a thousand tests. 83 positions differ significantly, out of the 583 that vary at
all.

## Slide 7: Which domains

So we know 83 positions differ. The question is whether they land anywhere
meaningful.

We grouped them by which part of the protein they fall in, and divided by how
long each part is, so a long region can't win just by being long.

Top of the list is the furin cleavage site. That's a five-residue stretch where
the protein has to be cut for the virus to fuse with a cell. Alpha changes it one
way, Delta changes it a different way. Two variants altering the same tiny region
by different routes suggests the same pressure acting twice.

Then the binding domain, then the N-terminal domain, which is a region antibodies
tend to target.

So the differences aren't random. They sit where the protein does its work.

## Slide 8: The ablation

Now the second question. Is the variant signal actually in the binding domain?

We trained a model on the full protein. Perfect accuracy, 1.000. Then we deleted
the binding domain entirely and ran it again. Still 1.000. Then we gave it only
the N-terminal domain. Still 1.000.

A perfect score sounds impressive but it isn't, for the reason I gave earlier.
The labels come from these mutations, so the model is reading its own answer key.
A baseline that always guesses the commonest variant gets 0.333, so the model is
doing something, but nothing surprising.

What is useful is the ablation. Removing the binding domain costs nothing. So we
measured it directly. We gave the model random positions and increased how many.

About fifty positions picked at random get you to 0.95. The information is
repeated all over the protein, not stored in one place.

## Slide 9: The confound

Now the geography question. And this is where we got it wrong first.

Our first model, predicting Europe against North America within Omicron, scored
0.94. We were pleased with that.

Then we checked which positions the model was actually using. The top ones were
N440K, Q498R, and a deletion at 69 to 70.

Those are not geographic markers. Those are the mutations that separate Omicron's
sub-types from each other, BA.1 from BA.2 from BA.5.

And here's why that helped it. BA.1 is about 31 per cent of our North American
samples, and under 3 per cent of our European ones. So a model that learns to
spot BA.1 gets most of the way to North America without using geography at all.

This is the same circularity I described at the start, one level down. We thought
we had removed it by holding the variant fixed. But a variant contains several
lineages, and those lineages aren't spread evenly across the world.

So we ran it again, holding the actual lineage fixed, not just the variant.

## Slide 10: The corrected result

And the signal survives.

Within B.1.1.7, which is Alpha, we get 0.659. Within Delta, 0.832. Within BA.2,
0.944. The baseline, which is just guessing the commoner region, is 0.500.

Three separate lineages, tested independently, all above baseline. So spike
sequence does carry information about where a sample came from.

Two honest caveats. We checked the collection dates for each lineage, so this
isn't just picking up that Europe and North America sampled at different times.

And there is a ceiling. Some identical sequences were collected on different
continents. No model that only sees the sequence can get all of those right. We
measured that limit at 0.871. So read the 0.659 against 0.871, not against a
perfect score.

## Slide 11: Limits

Four things we can't argue away.

First, and this is the one that bothers us most. Even inside a single lineage,
samples from the same place tend to be related to each other. So we can show the
signal is real and repeatable. We can't show that geography is what causes it.

Second, our data reflects where sequencing happened, not where the virus was.
After removing duplicates we had too few distinct sequences from Africa and Asia
to model them at all. So this result is Europe against North America, not the
world.

Third, our alignment method misplaces one deletion in Omicron, around positions
142 to 145. We've documented where that affects things.

And fourth, the variant prediction stays circular. We kept it as a check, not as
a result.

## Slide 12: Conclusions

So, three things.

The variants differ at a limited and consistent set of positions, 83 of them, and
they sit where the protein binds, where it gets cut, and where antibodies target
it.

That signal isn't localised. Deleting the binding domain costs nothing, and about
fifty arbitrary positions are enough to identify a variant.

And sequence carries partial information about geography. Above baseline inside a
fixed lineage, replicated three times, and bounded at 0.871.

If there's one thing I'd take from doing this, it's that twice a high score
turned out to be measuring the labels rather than the biology. Checking which
features the model was using is what caught it both times.

Everything is in the repository on the slide. Are there questions or feedback?

---

## During questions

Stay on slide 12. Don't advance.

Slide 14 has pre-built answers on the perfect score, leakage, and balanced
accuracy. Slide 13 is references, only if a source is queried.

If you don't know, say so and say where it is: "that's in Stage 5 of the
notebook, I can open it." That beats guessing.

**Don't say:** "100% accuracy" as a boast, always add that the task is circular.
"This proves geography causes the mutations," it doesn't. Any number not on a
slide or the cheat sheet.

**Expected question: why is the signal peptide not on your domain chart?**
It scores highest per residue, about 69 per 100. We left it off because it's a
13-residue stretch that gets cleaved off the mature protein, and it sits at the
start of the gene where sequences are often truncated. In our pipeline a
truncated start reads as deletions, so we think that number is a technical
artefact rather than biology. It's in the notebook with that caveat.
