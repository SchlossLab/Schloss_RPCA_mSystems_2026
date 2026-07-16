

\usepackage[left]{lineno}

\linenumbers

\modulolinenumbers

\usepackage{helvet}

\renewcommand*\familydefault{\sfdefault}

\usepackage[T1]{fontenc}

## Editor’s comments:

**The reviewers are in strong agreement that the current scope of the
study is too narrow and the presentation too one-sided for publication
in mSystems. Both recommend an expansion of the study to encompass both
more datasets and more metrics of success, and suggest adopting a more
neutral tone with a focus on when rarefaction is most appropriate and
when it may not be.**

**If these changes are beyond the scope of the study the manuscript may
be more appropriate for Spectrum but would still require some revision
and review.**

I appreciate the reviews provided by the two reviewers. Reading through
my responses to the reviewers you will find that I was not swayed by the
reviewers suggestions to make radical changes to my manuscript to
broaden its scope. Because Martino was originally published in
*mSystems*, I decided to submit my reproduction of their work to
*mSystems*. The goal of my study was to reassess Martino’s RPCA method
using their original simulation data and expand their analysis to
include rarefaction and other normalization and compositional data
analysis methods. This was made clear in the Abstract and Introduction
of the manuscript. Although the reviewers suggest a number of other
analyses and a desire for me to use more simulation frameworks, I
believe these are either well outside of the scope of what Martino did
or duplicative of two recent papers I published (Schloss 2023, 2024).
The reviewers sought a “comprehensive” study to settle the question of
rarefaction once and for all. This is an unrealistic expectation that
has never been achieve elsewhere in science.

I do not agree with the reviewers’ claims that my analysis was
one-sided. If anything I presented “other sides” by using the simulation
framework of proponents of compositional methods, replicated their
results, and then added a more comprehensive analytical framework by
using additional methods to controlling for uneven sampling effort. The
original Martino manuscript was so one-sided that I am not sure how they
could say I was being one-sided given my attempts to strengthen the
Martino study.

If you do not feel that my contributions are significant enough to be
published in the journal, I will understand. I ask that a decision be
rendered as soon as possible so that I can transfer the manuscript to
either *mSphere* or *Microbiology Spectrum*.

The line numbers I reference correspond to the marked up version of the
resubmitted manuscript.

## Reviewer \#1

**This manuscript revisits the long-standing debate on rarefaction vs
normalization/compositional approaches for uneven sequencing depth and
concludes that rarefaction is superior. The author re-analyzes Martino
et al.’s RPCA benchmarking, reports coding and methodological issues,
modifies the simulations, and re-runs comparisons, ultimately concluding
that rarefaction “remains the preferred method.” The author has made
valuable contributions in recent years advocating for rarefaction
through systematic reanalysis of prior benchmarking studies (Schloss,
2023, 2024.), establishing rarefaction’s practical utility. While the
paper raises important reproducibility concerns and provides useful
re-analysis, the current manuscript reads primarily as a targeted
rebuttal rather than a balanced and comprehensive benchmarking study.
Here are the major concerns**

As stated at L129-135, the goal of my paper was to use the benchmarking
framework originally used by Martino et al to compare their RPCA method
to other normalization and compositional analysis methods as well as
rarefaction. The reviewer is correct that this is a targeted rebuttal. I
contend that the analysis is balanced. I have reused Martino’s original
simulation framework and then improved the framework by fixing the
numerous coding errors present in their benchmarking code. This study
was not intended to be comprehensive. Rather it adds to the two papers
that the reviewer cited in this comment. A comprehensive benchmarking
study is well beyond the scope of a single study. It is more thorough
than Martino’s analysis, which cherry-picked analysis approaches that
are either not used in the microbiome field (e.g., KNN analysis) or that
put RPCA in a more favorable light (e.g., using datasets with large
effect sizes).

**1. The manuscript addresses a focused methodological question through
a reanalysis and modification of a previously published benchmark
(Martino et al.). While this re-evaluation is useful and technically
sound, the overall scope and impact may be limited. The study primarily
revisits an existing simulation framework, identifies implementation
issues, and evaluates rarefaction under revised settings using a
relatively narrow set of tasks (e.g., PERMANOVA power, KNN
classification, and sequencing-depth sensitivity). As currently framed,
the work reads more as a methodological clarification or refinement of a
prior benchmark than as a broadly impactful advance for the microbiome
field. To strengthen the manuscript’s impact, the evaluation should be
expanded to include a broader range of common microbiome analyses—such
as differential abundance testing, diversity estimation, distance
metrics, regression models, multivariate inference, and count-based
modeling—so that the conclusions are better supported across the wider
microbiome workflow.**

Because of the nature of RPCA, this really only makes sense in terms of
beta diversity analysis. The “relatively narrow set of tasks” are the
same tasks used in Martino et al. The other tasks the reviewer suggests
have been done by others and myself in previous studies. The current
study is already quite long to add other tools that do not necessarily
relate to downstream applications of RPCA. Furthermore, the focus of
this study was on the question of why I found one result using RPCA in
another simulation and benchmarking framework and Martino found a
different result.

**2. The manuscript is currently framed more as a rebuttal to a specific
prior study than as a neutral benchmark. The narrative repeatedly
focuses on identifying issues in Martino et al.’s work (e.g.,
implementation details, simulation choices, and comparison setup), which
shifts the focus toward a single-case critique and limits the breadth of
evidence needed to support general conclusions. A robust benchmarking
study would ideally include multiple methods, diverse simulation
frameworks, and evaluations across a range of assumptions to assess
generalizability. In contrast, the present work centers on critiquing
one benchmark, replacing it with a revised simulation framework, and
drawing broader conclusions. This narrower design limits the
generalizability of the findings and suggests that a more balanced,
multi-method benchmarking approach would strengthen the manuscript.**

I am not sure that there is such a thing as a “neutral benchmark”.
However, I used the same simulation and benchmarking strategy used by
Martino. I then improved it by fixing the numerous technical problems in
their implementation so that it more closely reflected the approach that
they described in their written text. In both cases, I found similar
results: rarefaction outperformed RPCA and other compositional methods.
These problems are not a matter of personal preference or style, they
are fundamental problems such as inadvertently transposing matrices (or
not), using nonsensical hyperparameters, and other issues that are
described in the detailed code review provided in the Supplemental Text
document. As outlined in my manuscript, this study is the third in a
series that have all used different benchmarking approaches. All have
come to a similar conclusion that rarefaction outperforms other methods
when applied to alpha and beta diversity metrics.

**3. The revised simulation framework still reflects a single modeling
perspective. Although the manuscript critiques the simulations of
Martino et al. and replaces them with a modified design, this does not
remove simulation dependence but rather substitutes one set of
assumptions for another. Key parameters, such as the number of features,
the number of differential taxa, and calibration of effect sizes to
achieve a specific PERMANOVA power, represent design choices that can
strongly influence comparative performance. A comprehensive benchmarking
study would ideally evaluate multiple simulation paradigms, including
established microbiome simulation tools (e.g., MIDASim, GUniFrac-based
simulations), along with a range of effect sizes, sparsity levels, and
compositional structures. As currently presented, the results primarily
reflect conclusions under a single simulation framework, which limits
the generalizability of the findings.**

My goal was not to remove simulation dependence in my re-analysis of the
RPCA study. The purpose of the revisions to the simulation framework
were to implement what the authors actually described in their
manuscript (e.g., homo and heteroscedastic noise) and to resolve clear
errors in their code. Using yet another framework to evaluate RPCA and
the other methods would not have answered the question of why Martino et
al found one result using RPCA and I found a different one. By using
Martino’s simulations with the same broad set of methods I used in
previous studies I was able to see that even with their original
simulation data I found that RPCA did not perform as well as other
methods. I was also able to identify coding errors in the code to
simulate and analyze their data.

**4. The manuscript demonstrates careful attention to detail in
examining the statement that “the most notable difference to Martino’s
version of the figure was that they inadvertently switched the labels of
the axes and performed the One CLR transformation on features rather
than samples.” However, this interpretation may be somewhat
over-scrutinized. Differences in data orientation are very common
between statistics and biology: statistical workflows typically place
samples in rows and variables in columns, whereas biological pipelines
often arrange features in rows and samples in columns. Depending on the
input format used in Martino’s implementation, their approach may simply
reflect the statistical convention and therefore may not require
transposing the data. This point should ideally be confirmed with the
authors of the original paper before drawing firm conclusions.**

This is not a matter of difference in data orientation. As stated in my
manuscript, the authors inadvertently switched the labels of the axes
and performed the One CLR transformation on features rather than
samples. This is elaborated upon in the Supplemental Text code review
where I point to the specific code on page 2 of the supplement. As
presented it appears that features overlap across samples. In reality,
samples overlapped across features - this doesn’t make sense.
Unfortunately, this error propagated through the rest of their
simulation data. By working through the code and data outputted at each
step, as I did in the Supplemental Text, it is clear this is not a
matter of statistical convention.

**5. This statement in line 155 appears overly critical. In practice,
focusing on the first two axes is a common and reasonable convention in
microbiome analyses because ordination plots are primarily used for
visualization. When group differences are strong, they typically appear
in the leading axes for PCA/PCoA-type methods. If signals emerge only in
higher dimensions, this more likely indicates a weaker or more diffuse
effect rather than a methodological flaw. RPCA does not inherently
require the use of two axes; the rank is a user-defined parameter
similar to the number of principal components in PCA. From this
perspective, the original choice to focus on two axes does not
constitute an error. Instead, the present manuscript’s extension to
additional axes can be viewed as a useful expansion that explores
additional scenarios. However, this extension is limited in scope and
largely represents a comparison with a single method, which reduces the
broader impact of this point.**

It is unclear to me what in L155 is overly critical. My sentence states
that my result agreed with their result, “Replicating what Martino
reported, the first two axes of the Bray-Curtis ordinations only
explained a fraction of the total variation in the data while the RPCA
ordinations explained all the variation.” As I go on to describe later
in the manuscript, increasing the number of axes for RPCA results in
significant changes to the analysis (e.g., Figure 7). I highlight that
there is no objective justification for 2 or 3 or 9 axes which may lead
to one selecting a number of axes that gives rise to a desired result
(L397-402).

**6. Some researchers view compositional transformations as an
alternative to rarefaction rather than a preprocessing step that must be
applied after rarefaction. In the compositional data analysis framework,
library size differences are treated as a scaling factor, and methods
such as CLR or other normalization approaches are intended to account
for this directly. From this perspective, applying rarefaction prior to
CLR may be considered redundant, as both steps are intended to address
sequencing depth effects, albeit from different methodological
viewpoints. This broader context may also help explain why Martino et
al. did not apply rarefaction prior to RPCA. Their work is grounded in
the compositional data analysis paradigm, which assumes that sequencing
counts are compositional and that differences in library size can be
addressed through normalization or log-ratio transformations rather than
subsampling. From this perspective, the absence of rarefaction is a
methodological choice consistent with that framework rather than an
oversight.**

I agree that combining rarefaction and compositional transformations is
redundant. I am unable to find where in the manuscript I criticized
Martino for not using rarefaction prior to RPCA. After all, The claim
made by Martino et al. is that those transformations replace rarefaction
and are not sensitive to differences in read depth. My criticism of
Martino was that for Bray-Curtis calculations they rarefied their data
to a 10-fold smaller depth than was necessary (L311-316). I went on to
point out that the difference didn’t really matter (L317-318).

**7. In line 345, the statement that RPCA-based methods performed worse
when the effect size was small should be interpreted cautiously. Reduced
performance in low–effect size settings may reflect the well-known
trade-off associated with log-ratio transformations, where variance
stabilization and compositional correction can come at the cost of
statistical power, particularly in sparse microbiome data. Reduced
performance in low–effect size settings is not unexpected and may
reflect the known trade-off associated with log-ratio transformations,
which tend to shrink variability and attenuate contrasts, particularly
in sparse microbiome data. As such, this observation does not
necessarily indicate a flaw specific to RPCA, but may instead reflect
broader properties of CLR-based approaches.**

Although these observations from the reviewer make sense to me, I am
unable to find a citation to them in the literature. If anything
McMurdie & Holmes argue that rarefaction causes a more significant hit
to power than the other methods they and I have used. If it is generally
true that rarefaction has more power than the other methods (as I and
the reviewer seem to agree), then I don’t see what the argument for
their use is since one does not know the effect size study *a priori*.
If one suspects they have a small effect size, then why use a
low-powered test? Because the analysis in this study is more focused on
RPCA than CLR-based methods in general, I will limit their critique in
this sentence to RPCA. I have added text at L393-395 that summarizes the
need for researchers to select the most powerful tests available.

**8. The rarefaction versus normalization debate has persisted for over
a decade, originating from fundamentally different statistical
philosophies. McMurdie & Holmes (2014) criticized rarefaction as
“inadmissible,” yet acknowledged that transformations can reduce power
depending on model assumptions—a critical nuance often ignored. Weiss et
al. (2017, Microbiome) demonstrated that normalization performance
depends strongly on data characteristics, with no single method
universally superior. Despite this established complexity, the current
manuscript fails to comprehensively address why this debate persists or
how to resolve it. The analysis does not answer fundamental questions:
Under what specific conditions (sparsity, effect size, library size
variation) does each approach maintain valid error control and optimal
power? How do rarefaction’s theoretical limitations (no variance model,
arbitrary thresholds, information loss) compare quantitatively to
log-ratio transformations’ known drawbacks in sparse data (Gloor et al.,
2017, Front. Microbiol.)? By selectively benchmarking RPCA against
rarefaction using metrics and datasets that favor depth-equalizing
methods, rather than systematically characterizing when different
approaches succeed or fail, the manuscript perpetuates advocacy rather
than providing principled guidance. Both approaches have documented
limitations. Rarefaction has been criticized for inducing bias and
discarding data (Willis, 2019; McMurdie & Holmes, 2014), while
normalization methods show data-dependent performance and can inflate
false discovery rates (Weiss et al., 2017; Lin & Peddada, 2020). A
comprehensive evaluation would focus on developing decision frameworks
for method selection based on dataset characteristics. The long-standing
debate will continue until the field moves beyond polarized positions
and establishes clear, evidence-based guidelines for when rarefaction’s
empirical robustness versus normalization’s statistical rigor is most
appropriate.**

The reviewer is correct that I don’t state “under what specific
conditions (sparsity, effect size, library size variation) does each
approach maintain valid error control and optimal power”. That’s beyond
the scope of this paper and I suspect beyond the scope of any single
study. My goal was to use the simulations implemented by Martino et
al. within a more comprehensive framework than they used to see why they
were so positive about their method and my previous analyses were less
so. Along the way, I found significant problems in their benchmarking
code. My analysis is far more comprehensive than the original Martino.
Taken with my two other papers that the reviewer had cited in an earlier
point, the analysis is now more comprehensive because it uses the
simulation framework proposed by Martino. It joins the conversation with
the other literature. I find the Reviewer’s comments to be unfair when
they appear to take the papers they cite as truth, when several of them
are reviews (e.g., Lin & Peddada, Gloor) or do not actually provide any
benchmarking (e.g., Willis). Meanwhile my previous papers appear to be
discounted in this discussion. Science moves forward one step at a time
correcting and adding to what has come before it. As an illustration of
this the version of the manuscript that I submitted cited all of the
studies that the reviewer cited. A single comprehensive study that
satisfies everyone has never been proposed in any area of science and is
unreasonable to expect here.

## Reviewer \#2

**This manuscript challenges the prevailing compositional-data paradigm
by arguing that rarefaction can outperform RPCA and related methods in
controlling uneven sequencing depth. The author presents extensive
benchmarking across the Martino et al. datasets and additional datasets,
which is valuable. However, several conceptual, methodological, and
presentation issues should be addressed before the conclusions can be
fully supported.**

**1. Although rarefaction outperforms compositional methods in certain
scenarios, the manuscript draws a very strong conclusion that
“Rarefaction remains the preferred method.” This framing is
unnecessarily polarizing and overlooks the contextual nature of
normalization choices. The more relevant question is not which method is
universally superior, but under which conditions each method is
appropriate. The author should explicitly discuss scenarios such as:**
**a) part of samples with low sequencing depth or sparse counts;** **b)
strong batch effects and cross-study heterogeneity in amplicon
datasets;** **c) variability in rRNA gene copy number across microbes.**

**A nuanced decision framework (rather than a universal recommendation)
would significantly improve the manuscript’s impact and fairness.**

The reviewer raises a number of interesting conditions that could be
included. Although I believe that their suggestion (a) is covered in
this and my previous analyses, I am reluctant and unsure of how to
account for (b) and (c). I am not familiar with arguments for using
compositional methods when there are batch effects or variation in rRNA
copy number. There would appear to be other approaches to correct the
underlying counts before downstream processing by rarefaction or any of
the other methods that I used. Because these conditions were not
included by Martino and my goal was to reproduce and strengthen their
analysis, I believe these considerations are outside of the scope of my
study. Frankly, in my three studies, I have yet to identify a nuanced
decision framework - rarefaction is the most reliable method for
controlling for uneven sampling effort. Other methods are at best as
good as rarefaction and often much worse. The safest choice would be to
use rarefaction regardless of the dataset.

**2. The manuscript primarily evaluates preprocessing methods using
classification and PERMANOVA. However, microbiome studies also aim to
identify biomarkers, infer ecological structure, or construct
co-occurrence networks. Rarefaction and CLR-based transformations can
affect these downstream analyses differently. At minimum, the author
should discuss potential implications for biomarker discovery and
network inference, or ideally provide benchmarking results for at least
one such downstream task.**

This study was designed to reproduce Martino and build upon it by adding
additional methods of correcting for multiple comparisons within an an
analysis framework based on beta diversity metrics. Other studies of
mine (Schloss 2023 and 2024) have looked at alpha and beta diversity as
well. Other studies that I cited look at biomarker discovery (e.g.,
Weiss et al). Considering the study is already quite long and biomarker
discovery or network inference were not included in Martino, I am
reluctant to add them to my reproduction of their analysis.

**3. In datasets with unequal group sizes (e.g., the sponge dataset),
accuracy alone is not an adequate performance metric. AUROC, AUPRC, or
balanced accuracy would provide a more comprehensive assessment. The
manuscript should report at least AUROC for imbalanced scenarios.**

In the analysis performed by Martino they used a balanced design by
randomly subsetting the number of samples from each treatment group.
Furthermore, the original Martino study used accuracy as the metric to
assess classification performance. To keep the analysis parallel between
my study and theirs, I think it is warranted to use accuracy. I don’t
think much is gained by using an alternative metric.

**4. The author states that weighted UniFrac was not included due to
difficulty in generating reliable trees and because Bray-Curtis and
weighted UniFrac rarely disagree (but no reference was cited in the
manuscript). Since all sequences were mapped to the SILVA database,
which provides a reference phylogeny. Therefore, UniFrac distances could
be computed using the SILVA tree.**

I have added the weighted UniFrac analysis requested by the reviewer
(introduced at L299-303). Through integrating this into my analysis, I
found an additional problem with Martino’s analysis of their case
studies. Namely, the trees that they supplied did not include all of the
sequences in their OTU table. As a result, they removed the missing
sequences from the OTU table. Unfortunately, the filtered OTU table was
not used to run RPCA. The result was another apples-to-oranges
comparison between the Bray-Curtis and weighted UniFrac analyses and the
RPCA analysis. Text describing this has been added to the Manuscript
(L283-285) and the Supplemental Text (P16, paragraph beginning with
“Third, in my attempts to reproduce the weighted UniFrac distances…”).

**5. The analysis focuses on amplicon OTU/ASV datasets. It remains
unclear whether the conclusions extend to shotgun metagenomic profiles,
which have different count distributions and compositional properties.
The author should either provide additional analyses or explicitly
discuss the expected applicability and limitations in shotgun
metagenomics.**

This manuscript only discusses amplicon sequence data. The only mention
of metagenomic data is in the opening sentence of the paper to state
that uneven sampling effort is a significant problem for both
approaches. I have opined on the need for rarefaction of metagenomics
data in other studies, but do not feel it is necessary to do so here
since this study is an attempt to reproduce and build upon the
simulations and analysis of Martino et al.

**6. For multi-panel figures (e.g., Figures 1 and 3), each panel should
be explicitly described in the figure legend. This will improve
readability and reduce ambiguity for readers and reviewers.**

Thank you for catching this oversight. The lettered panels have now been
described explicitly in the figure caption.

**7. This manuscript adopts a strongly adversarial tone toward Martino
et al., using phrases such as “straw man,” “critical errors,” and
“foreign to standard microbiome analyses.” While scientific critique is
appropriate, such language may be perceived as polemical and could
detract from the manuscript’s reception. I recommend moderating the
rhetoric and focusing on technical arguments rather than adversarial
framing.**

I suspect we all have different tolerances for tone of language. For
instance, “inadmissable” was used in the title of the McMurdie & Holmes
study that started the debate over rarefaction. I consider that word to
be needlessly inflammatory and incorrect. Regardless, I have removed
“straw man” from the text. I believe that “foreign to standard
microbiome analyses” was the most charitable way of describing analyses
that have yet to be used in microbiome studies (e.g., KNN
classification). I believe this is a more tolerable statement than
simply saying “wrong.” Similarly, the errors that I outlined in the
extensive detailed code review caused critical problems in the
simulations that caused them to deviate significantly from what the
authors described in the text of their papers. These errors were
critical.

**8. In addition, I also suggest to change the article type from
“Research Article” to “Methods and Protocols” or
“Opinions/Hypotheses”.**

I do not see how the manuscript I submitted would fall under either of
these categories. It is not a method or protocol and is not an opinion
or hypothesis. If *mSystems* has a format like *Nature*’s “Matters
Arising”, I would consider that but suspect it is too long for such a
format.
