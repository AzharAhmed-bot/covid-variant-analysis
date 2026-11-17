# Presentation Script

14 slides, six speakers, about 13 minutes. That leaves room to go slightly over
without running past a 15-minute slot.

Lines in **bold** are what to say. Text in brackets is a stage direction, not
something to read out. You do not have to say it word for word. Say it in your
own words, but keep the order and keep the numbers exact.

**Speaking order**

| Speaker | Slides | Topic | Time |
|---|---|---|---|
| Stephen Nyariki | 1–2 | Opening and why the question changed | 2:00 |
| Omar Siddique | 3–4 | What we asked and what we collected | 2:10 |
| Shadrack Njau | 5–6 | Mutation counts and where they fall | 2:10 |
| Hetvi Barot | 7–8 | Which domains, and the ablation | 2:20 |
| Casper Muhereza | 9–10 | The confound and the corrected result | 2:40 |
| Azhar Takoy | 11–12 | Limits, conclusions, questions | 2:00 |

---

## Stephen Nyariki

### Slide 1: Title (45 seconds)

**Good morning. We're group [X], and our project is on the SARS-CoV-2 spike
protein.**

*[Point at the red knobs on the image.]*

**Those red parts are the spike protein. It's what the virus uses to get into
human cells, and it's what changes most between variants.**

**Our finding is in the title. If you give us a virus sequence, telling you
which variant it is turns out to be easy. Telling you where in the world it
came from is much harder. The talk explains both halves of that.**

### Slide 2: Why the question changed (75 seconds)

**First, something we have to be upfront about. The question we got approved
is not the question we answered, and here's why.**

**Variants like Alpha and Delta aren't discovered by looking at the data. They
are defined by a list of mutations. A sequence is called Omicron because it
carries the Omicron mutations.**

*[Point at the table.]*

**So when our original question asked whether Omicron has more mutations in the
binding domain, the answer is right there in the catalogue. Alpha has one,
Delta has two, Omicron has fifteen to seventeen. We didn't need to collect
anything to know that.**

**And the second half, training a model to predict the variant, is the same
problem. The model would just be learning the rule that was used to apply the
labels in the first place.**

**So we kept that comparison, but only as a check that our pipeline works. The
real project is built on questions where we genuinely did not know the answer.**

**Omar will take those.**

---

## Omar Siddique

### Slide 3: The three questions (60 seconds)

**Thanks Stephen. We settled on three.**

**One: where along the protein do the variants actually differ? Not just at
the catalogue positions, but across all 1,273 of them.**

**Two: is that difference concentrated in the binding domain, which is the part
everyone focuses on? We tested that by deleting the binding domain from the data
and seeing how much accuracy we lost.**

**Three, and this is the one that matters: can a model tell which part of the
world a sample came from, using only the sequence?**

**That third one is the honest question, because where a sample was collected
has nothing to do with how a lineage is defined. So the model has to find real
signal. It could have found nothing, and we agreed beforehand that we'd report
that if it happened.**

### Slide 4: The data (70 seconds)

**On the data. We pulled 14,965 records from NCBI.**

*[Point down the numbers on the right.]*

**After we required a sequence and a usable collection date, and dropped
sequences with too many unresolved positions, we were left with 8,756.**

**Then the big cut. 1,406. That's because most samples of the same variant are
identical in spike. Thousands of records collapse to a few hundred genuinely
different sequences.**

**We had to remove those duplicates, because if the same sequence appears in
both the training and the test set, the model can score well just by
remembering it.**

**And that was the real limit on this project. Not how much we could download,
but how many distinct sequences exist. We went from 2,600 records to 15,000 and
roughly doubled the distinct sequences, not six times.**

**Shadrack has the first results.**

---

## Shadrack Njau

### Slide 5: Mutation counts (70 seconds)

**This is the check Stephen mentioned.**

*[Point at the two panels.]*

**On the left, total mutations per sample. On the right, only those in the
binding domain. Each dot is one virus sample.**

**The medians are one for Alpha, two for Delta, fifteen for Omicron.**

**Now, those are the published numbers. But we never gave the catalogue to our
code. It worked them out from the sequences alone and landed on the same
values.**

**That matters, because if this slide had come out wrong, nothing after it
would be worth showing you. The statistics are at the bottom: the difference
between the three groups is significant, and Omicron separates completely from
the other two.**

### Slide 6: Where the changes fall (60 seconds)

**Now where along the protein those changes sit. The axis runs across all
1,273 positions, one row per variant.**

*[Point at the Omicron row.]*

**Two things to notice. The changes aren't spread out evenly. Each variant has
a small set of positions that nearly every sample has changed, and long
stretches that never change at all.**

**And Omicron has that dense band through the shaded region, which is the
binding domain. Alpha and Delta don't.**

**We tested every position separately, and corrected for the fact that we ran
over a thousand tests. 83 positions differ significantly out of the 583 that
vary at all.**

**Hetvi will take it from here.**

---

## Hetvi Barot

### Slide 7: Which domains (70 seconds)

**Thanks. So we know 83 positions differ. The question is whether they land
anywhere meaningful.**

**We grouped them by which part of the protein they fall in, and divided by how
long each part is, so a long region can't win just by being long.**

*[Point at the top bar.]*

**Top of the list is the furin cleavage site. That's a five-residue stretch
where the protein has to be cut for the virus to fuse with a cell. Alpha changes
it one way, Delta changes it a different way. Two variants altering the same
tiny region by different routes suggests the same pressure acting twice.**

**Then the binding domain, then the N-terminal domain, which is a region
antibodies tend to target.**

**So the differences aren't random. They sit where the protein does its work.**

### Slide 8: The ablation (70 seconds)

**Now the second question. Is the variant signal actually in the binding
domain?**

*[Point at the table.]*

**We trained a model on the full protein. Perfect accuracy, 1.000. Then we
deleted the binding domain entirely and ran it again. Still 1.000. Then we gave
it only the N-terminal domain. Still 1.000.**

**Now, a perfect score sounds impressive but it isn't, for the reason Stephen
gave. The labels come from these mutations, so the model is reading its own
answer key. A baseline that always guesses the commonest variant gets 0.333, so
the model is doing something, but nothing surprising.**

**What is useful is the ablation. Removing the binding domain costs nothing. So
we measured it directly: we gave the model random positions and increased how
many.**

*[Point at the curve.]*

**About fifty positions picked at random get you to 0.95. The information is
repeated all over the protein, not stored in one place.**

**Casper has the part we're most pleased with.**

---

## Casper Muhereza

### Slide 9: The confound (90 seconds, take your time)

**Thanks. So, the geography question. And this is where we got it wrong first.**

**Our first model, predicting Europe against North America within Omicron,
scored 0.94. We were pleased with that.**

**Then we checked which positions the model was actually using. The top ones
were N440K, Q498R and a deletion at 69 to 70.**

**Those are not geographic markers. Those are the mutations that separate
Omicron's sub-types from each other, BA.1 from BA.2 from BA.5.**

*[Point at the BA.1 row.]*

**And here's why that helped it. BA.1 is about 31 per cent of our North
American samples, and under 3 per cent of our European ones. So a model that
learns to spot BA.1 gets most of the way to "North America" without using
geography at all.**

**This is the same circularity Stephen described at the start, one level down.
We thought we'd removed it by holding the variant fixed. But a variant contains
several lineages, and those lineages aren't spread evenly across the world.**

**So we ran it again, holding the actual Pango lineage fixed, not just the
variant.**

### Slide 10: The corrected result (70 seconds)

**And the signal survives.**

*[Point at the three numbers.]*

**Within B.1.1.7, which is Alpha, we get 0.659. Within Delta, 0.832. Within
BA.2, 0.944. The baseline, which is just guessing the commoner region, is
0.500.**

**Three separate lineages, tested independently, all above baseline. So spike
sequence does carry information about where a sample came from.**

**Two honest caveats on this slide. We checked the collection dates for each
lineage, so this isn't just picking up that Europe and North America sampled at
different times.**

**And there's a ceiling. Some identical sequences were collected on different
continents. No model that only sees the sequence can get all of those right. We
measured that limit at 0.871. So read the 0.659 against 0.871, not against a
perfect score.**

**Azhar will close.**

---

## Azhar Takoy

### Slide 11: Limits (60 seconds)

**Thanks Casper. Four things we can't argue away.**

**First, and this is the one that bothers us most. Even inside a single
lineage, samples from the same place tend to be related to each other. So we
can show the signal is real and repeatable. We can't show that geography is
what causes it.**

**Second, our data reflects where sequencing happened, not where the virus was.
After removing duplicates we had too few distinct sequences from Africa and
Asia to model them at all. So this result is Europe against North America, not
the world.**

**Third, our alignment method misplaces one deletion in Omicron, around
positions 142 to 145. We've documented where that affects things.**

**And fourth, the variant prediction stays circular. We kept it as a check, not
as a result.**

### Slide 12: Conclusions (60 seconds)

**So, three things.**

**The variants differ at a limited and consistent set of positions, 83 of them,
and they sit where the protein binds, gets cut, and gets targeted.**

**That signal isn't localised. Deleting the binding domain costs nothing, and
about fifty arbitrary positions are enough to identify a variant.**

**And sequence carries partial information about geography. Above baseline
inside a fixed lineage, replicated three times, and bounded at 0.871.**

**If there's one thing we'd take from doing this, it's that twice a high score
turned out to be measuring the labels rather than the biology. Checking which
features the model was using is what caught it both times.**

**Everything is in the repository on the slide. Are there questions or
feedback?**

*[Stay on this slide. Do not advance.]*

---

## During questions

- Azhar takes questions and passes them to whoever owns that slide.
- Slide 14 is a pre-built answer slide. Go to it if asked about the perfect
  score, leakage, or balanced accuracy.
- Slide 13 is references. Only if a source is queried.
- If nobody knows, say so and say where it is: "that's in Stage 5 of the
  notebook, I can open it." That beats guessing.
- The cheat sheet has twelve likely questions with short answers. Have it open.

## Do not say

- "100% accuracy" as a boast. Always add that the task is circular.
- "This proves geography causes the mutations." It doesn't.
- Any number that isn't on a slide or the cheat sheet.

## One question we expect

**Why is the signal peptide not on your domain chart?**

*Answer:* It actually scores highest per residue, about 69 per 100. We left it
off because it's a 13-residue stretch that gets cleaved off the mature protein,
and it sits at the start of the gene where sequences are often truncated. In our
pipeline a truncated start reads as deletions, so we think that number is a
technical artefact rather than biology. It's in the notebook with that caveat.
