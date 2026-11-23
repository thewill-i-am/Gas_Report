# Grover Adaptive Search for Biomarker Selection

## 1. What are we searching for?

We search for a **small set of biomarkers**, not for a participant or an autism diagnosis.

For each dataset, let:

- \(N\): number of candidate biomarkers.
- \(k\): number of biomarkers to select, here \(k\in\{3,5,7\}\).
- \(x_i\in\{0,1\}\): whether candidate biomarker \(i\) is selected.
- \(x=(x_1,\ldots,x_N)\): a candidate panel, with \(\sum_i x_i=k\).

| Dataset | Candidates \(N\) | Meaning of a selected feature |
|---|---:|---|
| S1 | 24 | A biochemical measurement |
| Plasma | 56 | A plasma metabolite/lipid measurement |
| Shotgun | 54 | A microbial abundance feature |

**Pipeline**

```text
Measurements and known ASD/TD labels
    → classical relevance and redundancy statistics
    → QUBO coefficients
    → GAS selects a panel of exactly k biomarkers
    → LDA / SVM uses those measurements to predict ASD or TD
    → cross-validated predictions give balanced accuracy
```

The statistical information enters the search through the QUBO coefficients and,
for guided GAS, the preparation weights. We do not amplitude-encode each
participant's measurement vector in this search.

**Implementation status:** the report uses ideal analytical amplitude-amplification
sampling, not execution of the full reversible oracle on quantum hardware.
The registers and circuits below describe a coherent realization of that search.
They are not a claim that these gates were executed for the reported benchmarks.

## 2. The QUBO: what makes a panel good?

QUBO means **Quadratic Unconstrained Binary Optimization**. Our reported GAS
variants use a quadratic objective but restrict the search to fixed-size panels:

$$
\min_{\substack{x\in\{0,1\}^N\\\sum_i x_i=k}} E(x),
\qquad E(x)=x^\top Qx.
$$

Let \(I_i\) be normalized relevance and \(R_{ij}\) pairwise redundancy, with
\(R_{ii}=0\). For the mean-normalized objectives in this report:

$$
Q_{ii}=-\frac{\alpha I_i}{k},
\qquad
Q_{ij}=\frac{(1-\alpha)R_{ij}}{k(k-1)}
\quad(i\ne j).
$$

For a selected set \(S\), this is equivalent to maximizing

$$
F(S)=-E(S)
=\frac{\alpha}{k}\sum_{i\in S}I_i
-\frac{1-\alpha}{\binom{k}{2}}
\sum_{\substack{i<j\\i,j\in S}}R_{ij}.
$$

In plain English: **reward informative features and penalize redundant pairs**.

### Avoid double-counting pairs

For a symmetric matrix:

$$
x^\top Qx=\sum_i Q_{ii}x_i+2\sum_{i<j}Q_{ij}x_ix_j.
$$

Thus a canonical coefficient stored once per pair is \(2Q_{ij}\), not \(Q_{ij}\).

### Where do the coefficients come from?

For `pearson_mean`, used by the selected Plasma \(k=3\) configurations:

$$
r_i=|\operatorname{corr}(X_i,y)|,\qquad
I_i=\frac{r_i}{\max_j r_j},\qquad
R_{ij}=|\operatorname{corr}(X_i,X_j)|,\quad i\ne j.
$$

Here \(y\) is the known ASD/TD label. Zero maxima are handled with positive
denominator floors in the code. For this objective, \(\alpha=0.77\).

Other tested objectives change the relevance statistic, the redundancy statistic,
or its weight: Fisher, mutual information, Spearman, LDA contribution, and ReliefF.
`fisher_weak` uses \(\alpha=0.95\); `fisher_only` uses \(\alpha=1\) and no redundancy.
These are different objectives, not merely different random seeds.

### Fixed size is enforced by the search space

The reported GAS search contains only panels with exactly \(k\) selected features.
It does not need a cardinality penalty in its objective.

An alternative search over all \(2^N\) bitstrings could add
\(\lambda(\sum_i x_i-k)^2\) or explicitly reject infeasible strings. That is not
the fixed-size preparation used in this report.

### Integer arithmetic for a coherent oracle

A reversible circuit needs a specified finite-precision representation.
The code quantizes canonical coefficients at scale \(L=2^{20}\):

$$
c_i=\operatorname{round}(LQ_{ii}),\qquad
c_{ij}=\operatorname{round}(2LQ_{ij})\quad(i<j).
$$

The integer energy is

$$
E_{\mathrm{int}}(x)=\sum_i c_ix_i+\sum_{i<j}c_{ij}x_ix_j.
$$

The search implementation uses the opposite sign:

$$
F_{\mathrm{int}}(x)=-E_{\mathrm{int}}(x).
$$

**Lower energy and higher score mean the same thing.** Comparisons use the
quantized coefficients, not a newly rounded floating-point score per panel.
Quantization can merge nearly equal objective values.

## 3. What is a register?

A **quantum register** is a group of qubits used together to represent a quantity.
A computational-basis state represents a bit pattern; a superposition can assign
amplitudes to many such patterns.

### Logical registers for a compute–compare–uncompute oracle

| Register | Size | Purpose |
|---|---|---|
| Selection register \(X\) | \(N\) qubits in a direct feature-mask encoding | Stores which biomarkers are selected |
| Score register \(V\) | \(b\) qubits | Temporarily stores the signed integer score |
| Comparison flag \(C\) | 1 qubit in this construction | Stores whether the score beats the threshold |
| Workspace \(W\) | Implementation-dependent | Temporary carries, pair products and comparator work |

The threshold \(T\) is a **classical value**, updated between rounds. It can be
compiled into arithmetic/comparator gates; it need not occupy a separate quantum
register. A design that loads \(T\) into a quantum register would require extra
qubits.

The initial logical state is

$$
|\psi\rangle_X|0\rangle_V^{\otimes b}|0\rangle_C|0\rangle_W.
$$

#### Selection-register example

For four candidate biomarkers \(A,B,C,D\), use the written order
\(|x_Ax_Bx_Cx_D\rangle\). Then:

$$
|1010\rangle_X\quad\longleftrightarrow\quad\{A,C\}.
$$

The reported implementation labels feature \(i\) with integer bit \(i\):
\(\mathrm{index}(x)=\sum_{i=0}^{N-1}2^ix_i\). This little-endian integer convention
should not be confused with the left-to-right feature order used in this example.

For Plasma, a direct mask register would contain **56 selection qubits**.
Selecting only three biomarkers does **not** mean using only three selection qubits.
There are \(\binom{56}{3}=27{,}720\) feasible panels in that register's fixed-weight
subspace.

A compact panel-index representation could use
\(\lceil\log_2\binom Nk\rceil\) qubits, but then the oracle needs reversible
decoding/unranking logic. It is not interchangeable with the direct mask encoding
without that extra construction.

#### How many score qubits?

Choose \(b\) so all intermediate signed arithmetic fits in

$$
[-2^{b-1},\,2^{b-1}-1].
$$

A conservative coefficient bound is

$$
B=\sum_i|c_i|+\sum_{i<j}|c_{ij}|.
$$

If explicitly computing \(F_{\mathrm{int}}-T\), a sufficient conservative choice
satisfies \(2^{b-1}>B+|T|\). Tighter bounds can reduce the width.
The scale \(2^{20}\) does **not** imply a 20-qubit score register: summation,
sign and subtraction range also matter. Overflow would invalidate the oracle.

## 4. State preparation: which panels can be sampled?

Let \(\Omega_k=\{x:\sum_i x_i=k\}\).

### Uniform GAS

$$
|\psi_u\rangle
=\frac{1}{\sqrt{\binom Nk}}\sum_{x\in\Omega_k}|x\rangle.
$$

Each feasible panel starts with probability \(1/\binom Nk\). This is a
fixed-weight (Dicke-type) state, not simply Hadamard gates on every qubit, which
would also include panels of the wrong size.

### Guided GAS

The implementation computes logits

$$
\ell_i=\operatorname{clip}\left(
4\left[I_i-\frac{1}{N-1}\sum_{j\ne i}\widehat R_{ij}\right],
-12,12\right),
\qquad
\widehat R=\frac{R}{\max_{a,b}R_{ab}}.
$$

When the redundancy maximum is zero, the normalized redundancy is zero.
The fixed-size preparation probability is

$$
p_g(S)=
\frac{\exp(\sum_{i\in S}\ell_i)}
{\sum_{T:\,|T|=k}\exp(\sum_{j\in T}\ell_j)},
\qquad
|\psi_g\rangle=\sum_{S:\,|S|=k}\sqrt{p_g(S)}\,|S\rangle.
$$

This favors relevant, less redundant features before amplification.
It still has support on every feasible panel. It is a product-weight distribution
**conditioned on size \(k\)**, not independent single-qubit rotations without a
cardinality constraint.

The guide changes the starting distribution. For a fixed QUBO and threshold,
it does not change which panels the oracle marks. The fixed guide strength is 4;
changing \(\alpha\) does not directly change that strength.

## 5. What is the oracle?

An **oracle** is a reversible operation that tests a candidate against a specified
condition. Here the condition is:

$$
g_T(x)=
\begin{cases}
1,&F_{\mathrm{int}}(x)>T,\\
0,&F_{\mathrm{int}}(x)\le T.
\end{cases}
$$

The condition is strict: equal-score panels are not improvements.
The phase oracle acts on feasible states as

$$
O_T|x\rangle=(-1)^{g_T(x)}|x\rangle.
$$

It flips the phase of improving panels. **It does not directly increase their
measurement probability.** Probability changes after the reflection/amplification
step.

### Reversible implementation

1. **Compute the score**
   \[
   |x\rangle|0\rangle_V\longrightarrow|x\rangle|F_{\mathrm{int}}(x)\rangle_V.
   \]
   Linear terms use additions controlled by \(x_i\); pair terms use additions
   controlled by \(x_ix_j\).
2. **Compare with \(T\)** and write the result to flag \(C\).
3. **Apply \(Z\) to the flag**, producing phase \((-1)^{g_T(x)}\).
4. **Uncompute the comparison** to return the flag and its workspace to zero.
5. **Uncompute the score** to return \(V\) and its workspace to zero.

Net action:

$$
|x\rangle|0\rangle_V|0\rangle_C|0\rangle_W
\longmapsto
(-1)^{g_T(x)}|x\rangle|0\rangle_V|0\rangle_C|0\rangle_W.
$$

Uncomputation removes information about \(x\) from the ancillary registers so the
selection register can interfere correctly. A phase-kickback design with an
ancilla in \(|-\rangle\) is another implementation; it is not an additional
mandatory register on top of this flag-and-\(Z\) construction.

**The oracle does not train LDA/SVM, compute balanced accuracy, or diagnose ASD.**
It compares the fixed QUBO score against a threshold.

## 6. Grover amplification and adaptive search

If \(A|0\rangle=|\psi\rangle\), the reflection about the prepared state is

$$
D_\psi=2|\psi\rangle\langle\psi|-I
=A(2|0\rangle\langle0|-I)A^\dagger.
$$

One Grover iteration is \(G_T=D_\psi O_T\).
For guided preparation the reflection must use that guided state, not a uniform
diffuser.

If the initial probability of improving panels is
\(a=\sum_{x:g_T(x)=1}p(x)\) and \(\theta=\arcsin\sqrt a\), then after \(j\) iterations:

$$
P(\text{improvement})=\sin^2((2j+1)\theta).
$$

Amplification is not monotonic in \(j\): too many iterations can reduce success.

### Policy used in the report

- Begin with a uniformly sampled feasible incumbent panel.
- Set \(T\) to its integer score.
- Perform 32 rounds, with a maximum of one Grover iteration per round.
- The search window starts at 1, so the first round uses \(j=0\).
- After failure the window becomes 2, allowing random \(j\in\{0,1\}\).
- After improvement, replace the incumbent and threshold, and reset the window to 1.
- Uniform or guided preparation supplies the distribution for each round.

```text
best = uniformly sampled size-k panel
T = integer_score(best)
window = 1
repeat 32 times:
    j = uniformly sampled integer from {0, ..., min(window, 2)-1}
    candidate = sample the ideal distribution after j Grover iterations
    if integer_score(candidate) > T:
        best = candidate
        T = integer_score(candidate)
        window = 1
    else:
        window = min(2 * window, 2)
return best
```

The benchmark simulates these probabilities using enumerated/ranked panel scores
and CDF sampling. It does not obtain a hardware speedup by doing so, and it does
not guarantee finding the globally optimal panel.

## 7. Worked example: four candidates, choose two

This is an **illustrative toy problem**, not measured ASD data.

Take \(N=4\), \(k=2\), \(\alpha=0.8\), relevance

$$
I=(1.0,\ 0.8,\ 0.6,\ 0.2),
$$

and redundancy

$$
R=\begin{pmatrix}
0&0.9&0.1&0\\
0.9&0&0.1&0\\
0.1&0.1&0&0.1\\
0&0&0.1&0
\end{pmatrix}.
$$

The symmetric QUBO matrix is

$$
Q=\begin{pmatrix}
-0.40&0.09&0.01&0\\
0.09&-0.32&0.01&0\\
0.01&0.01&-0.24&0.01\\
0&0&0.01&-0.08
\end{pmatrix}.
$$

For two selected features:
\(F(\{i,j\})=0.4(I_i+I_j)-0.2R_{ij}\).

| Panel | Mask \((x_A,x_B,x_C,x_D)\) | Score \(F=-E\) | Integer score, toy scale 100 |
|---|---|---:|---:|
| A, B | 1100 | 0.54 | 54 |
| A, C | 1010 | **0.62** | **62** |
| A, D | 1001 | 0.48 | 48 |
| B, C | 0110 | 0.54 | 54 |
| B, D | 0101 | 0.40 | 40 |
| C, D | 0011 | 0.30 | 30 |

Although B is individually more relevant than C, A and B are highly redundant.
The objective therefore prefers **A+C** to A+B.

For transparent arithmetic use scale 100 in this toy example, not the production
scale \(2^{20}\). The score is

$$
F_{\mathrm{int}}=
40x_A+32x_B+24x_C+8x_D
-18x_Ax_B-2x_Ax_C-2x_Bx_C-2x_Cx_D.
$$

At demonstration threshold \(T=55\), only A+C is marked:

$$
O_{55}|1010\rangle=-|1010\rangle.
$$

All other feasible basis states retain their sign. A signed 8-bit accumulator
is sufficient for this toy score and subtraction by 55 using its coefficient
bounds; ancillary workspace is additional.

With uniform preparation, the initial success probability is \(1/6\). After one
Grover iteration:

$$
\sin^2\left(3\arcsin\sqrt{\frac16}\right)
=\frac{49}{54}\approx90.74\%.
$$

**90.74% is the probability of sampling the marked panel in this toy search.
It is not an ASD classification accuracy.**

The threshold 55 illustrates the oracle. In the actual adaptive policy the
threshold comes from the incumbent: starting with A+B gives \(T=54\), which
also marks only A+C because the comparison is strict. If A+C is sampled, the
threshold becomes 62; no feasible panel in this toy problem improves on it.

## 8. How selection relates to ASD classification

Once GAS returns a panel, the measured values of those biomarkers are inputs to
LDA or linear SVM. Classifier LOOCV trains on all but one participant and predicts
the held-out participant's ASD/TD label, repeated across the cohort.

$$
\text{Balanced accuracy}
=\frac12\left(
\frac{TP}{TP+FN}+\frac{TN}{TN+FP}
\right).
$$

ASD is the positive class: the first term is sensitivity, and the second is
specificity for TD.

Keep three quantities separate:

| Quantity | Role |
|---|---|
| QUBO score | Determines which panels GAS marks as improvements |
| Sampling probability | Determines how likely GAS is to return a candidate |
| Balanced accuracy | Evaluates ASD/TD classification using the selected measurements |

We rank tested objective configurations by mean LDA/SVM balanced accuracy over
20 seeds. For the presentation's “top panel,” we then select the highest-scoring
observed seed panel within the winning configuration, using mean LDA/SVM BA.
That selected maximum is distinct from the configuration average.

The QUBO and guide in these experiments use the full cohort before classifier
LOOCV, and objective selection is retrospective. The results are exploratory,
not independent validation of an autism diagnostic test.

## 9. Thirty-second explanation

> “Each search state represents a combination of biomarkers. The QUBO rewards
> relevance to ASD versus TD and penalizes redundant features. A threshold oracle
> marks combinations better than the current panel, and Grover amplification
> increases their sampling probability. We then use LDA and SVM to evaluate how
> well the selected measurements classify participants. Balanced accuracy is the
> evaluation and objective-selection criterion, not the score inside the oracle.”

## 10. Implementation references

Paths are relative to the repository root:

- `poc/experiments/plasma_gas_objectives/gas_objectives.py`: mean-normalized objectives and the report's GAS campaign.
- `poc/experiments/GAS_vs_QAOA/objective.py`: coefficient quantization, statistics, guide logits and fixed-size encoding.
- `poc/experiments/plasma56_benchmark/run_plasma.py`: signed integer score representation and 32-round cap-1 policy.
- `poc/experiments/biomarker_count_benchmark/`: fixed-cardinality benchmark and analytical sampling.
- `presentations/gas_objective_round2/build_top1.py`: objective ranking, report descriptions and best-observed panels.
- `poc/TECNICA_GROVER_QUBO.md`: broader reversible-oracle design documentation; distinct from claiming hardware execution of this report.
