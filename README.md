# Minimum-weight codewords of Polar, PAC, PS-PAC and CRC-polar codes

This directory contains two versions of the same Python program:

- `mwc_enumerator.py`: Python standard library only.
- `mwc_enumerator_numba.py`: the same calculation, with Numba used for the time-consuming loops.

The program finds the minimum distance $d_{\min}$ and the number of minimum-weight codewords $A_{d_{\min}}$ for **one code at a time**. Set `code_type` to `polar`, `pac`, `ps-pac` or `crc-polar`. The minimum distance is an output, not a parameter that must be known beforehand.

The usual calculation is based on the structure of minimum-weight codewords in polar and PAC codes. If this calculation does not determine $d_{\min}$, the program can also count low-weight words using a parity-check matrix. That search is described separately below. If it cannot be completed within its limits, the program does not assign an unverified value to $A_{d_{\min}}$.

## Running the program

Python 3.10 or later is required. The first script has no additional dependencies. To use the Numba version, install NumPy and Numba:

```bash
python -m pip install numpy numba
```

The simplest way to specify a code is to edit `CONFIG` near the beginning of the script. For example, for a $(64,32)$ PS-PAC code:

```python
CONFIG = dict(
    code_type="ps-pac",
    N=64,
    K=32,
    construction="dega",
    design_snr_db=2.0,
    alpha=7,
    pac_poly="1011011111",
)
```

Run **either** version:

```bash
python mwc_enumerator.py
```

or

```bash
python mwc_enumerator_numba.py
```

The same parameters can instead be supplied on the command line. For example:

```bash
python mwc_enumerator.py --code-type ps-pac --N 64 --K 32 --construction dega --design-snr-db 2 --alpha 7 --pac-poly 1011011111
```

For the accelerated version, change the filename to `mwc_enumerator_numba.py`. When command-line arguments are supplied, they are used in place of `CONFIG` (they do not inherit values from the dictionary). Run either script with `--help` to see the available arguments. The `enumerate_code(...)` function can also be called from another Python file.

## Code parameters

For all four families, $N$ is the block length and $K$ is the number of independent information bits. For CRC-polar codes, **$K$ does not include the CRC bits**. A CRC polynomial of degree $t$ gives $K+t$ selected polar input positions.

`construction` specifies how the information set $\mathcal I$ is selected:

| Value | Construction |
| --- | --- |
| `dega` | Gaussian approximation (GA), with `design_snr_db` and, if needed, `design_rate` [7] |
| `crossed-rm` | Crossed RM-polar profile: all rows above the RM weight boundary, with the remaining positions chosen by GA reliability [6] |
| `5g` | Reliability sequence in 3GPP TS 38.212 [8] |
| `pw` | Polarization-weight ordering, with optional `pw_beta` [9, 10] |
| `explicit` | A user-specified information set, given by `info_set` |

The `5g` option uses the standard reliability sequence for $N\le 1024$; it does not carry out 5G rate matching or interleaving. All indices start at zero and refer to the rows of $G_N=F^{\otimes n}$, where $F=\left[\begin{smallmatrix}1&0\\1&1\end{smallmatrix}\right]$. No extra bit-reversal permutation is applied.

For an explicit information set, write indices and ranges as, for example, `"15,23,26-31,38,39"`. For `ps-pac`, an explicitly supplied information set is treated as *already shifted*, unless a full `reliability_order` is also given. In the latter case, the program applies the shift to the supplied base set.

Other parameters are:

| Parameter | Purpose |
| --- | --- |
| `pac_poly` | PAC convolutional coefficients, in increasing delay order; for example, `1011011` means $(1,0,1,1,0,1,1)$ |
| `alpha` | Number of positions exchanged in PS-PAC construction |
| `crc_poly` | CRC polynomial, such as `"x^5+x^3+1"` or `"0x29"` |
| `crc_positions` | CRC input indices; by default, the last $t$ selected positions |
| `per_coset` | Print the individual leading-coset contributions |

The convolutional precoder is defined by

$$u_i=\sum_{\ell=0}^{s}p_\ell v_{i-\ell}\pmod 2, \qquad v_j=0\quad(j<0).
$$

Here $v$ is the rate-profile vector and $u$ is the polar-transform input. For CRC-polar codes, the CRC remainder is calculated as $r(x)=d(x)x^t\bmod q(x)$, with the remainder inserted at the selected CRC positions in increasing index order. The CRC positions must occur after the data positions.

## Enumeration methods

### 1. Polar codes

The number of minimum-weight codewords is calculated using the row-based construction of Rowshan, Dau and Viterbo [1]. For a polar information set $\mathcal I$ satisfying the usual partial order, the weight of generator row $g_i$ is

$$w(g_i)=2^{\operatorname{wt}(i)},$$

where $\operatorname{wt}(i)$ is the number of ones in the binary expansion of $i$. Define

$$d=\min_{i\in\mathcal I}w(g_i), \qquad
\mathcal B=\{i\in\mathcal I:w(g_i)=d\},$$

and

$$\mathcal K_i=\{j\in\mathcal I:j>i,\ |\operatorname{supp}(j)\setminus\operatorname{supp}(i)|=1\}.$$

Then

$$d_{\min}=d, \qquad A_{d_{\min}}=\sum_{i\in\mathcal B}2^{|\mathcal K_i|}.$$

The program forms $\mathcal B$ and $\mathcal K_i$ directly from the information set and evaluates this sum. This is the Rowshan–Dau–Viterbo formula [1]; its minimum-weight count agrees with the earlier result of Bardet, Drăgoi, Otmani and Tillich for decreasing monomial codes [2]. If the information set fails the assumed partial-order condition, the program does not apply the formula as an exact multiplicity count.

### 2. PAC codes

The PAC calculation follows the minimum-weight codeword construction and the capable-coset analysis in [1, 3, 4]. A polar minimum-weight codeword with leading row $g_i$ is described by

$$\mathcal O=\{i\}\cup\mathcal J\cup\mathcal M(\mathcal J),
\qquad \mathcal J\subseteq\mathcal K_i,$$

where $\mathcal M(\mathcal J)$ contains the balancing-row indices required to keep the weight unchanged. These indices are updated when a new element is added to $\mathcal J$, using the `addToM` rule of [3, 4]. Repeated contributions cancel over $\mathbb F_2$.

For each leading coset, the program examines frozen indices

$$\mathcal F_i^*=\{f\in\mathcal I^c:f>i,\ |\operatorname{supp}(f)\setminus\operatorname{supp}(i)|>1\}.$$

When $\mathcal F_i^*$ is empty, the coset is *incapable* of reducing the polar minimum-weight count, as described in [4]. Otherwise, put $f_{\max}=\max\mathcal F_i^*$ and enumerate the choices in the truncated core

$$\mathcal K_i^{\mathrm{tr}}=\{j\in\mathcal K_i:j<f_{\max}\}.$$

For each choice, the routine proceeds through the input positions in increasing order. It uses the convolution relation to determine the required value of $v_i$ at information positions; at frozen positions, $v_i=0$, so the resulting $u_i$ must agree with the polar MWC row combination. A choice that violates a frozen-position condition is rejected. The remaining core positions contribute a factor $2^{|\mathcal K_i|-|\mathcal K_i^{\mathrm{tr}}|}$, as in the truncated-core treatment [3, 4]. Counts are summed over the leading cosets.

This is a coset enumeration, not a search through all $2^K$ PAC messages. The `max_structural_subsets` parameter (default $2^{24}$) limits the number of truncated-core choices considered for a coset. It is a limit on **this structural calculation**, not on the syndrome search discussed below.

The resulting structural count applies to the minimum polar-row-weight class. Counts associated with higher row-weight leaders are not, by themselves, a complete higher-weight spectrum. If the lowest class is absent, an additional minimum-distance calculation is required.

### 3. Profile-shifted PAC (PS-PAC) codes

PS-PAC codes use the rate-profile modification proposed by Gu, Rowshan and Yuan [5]. First, a base information set of size $K$ is constructed. The `alpha` most reliable positions in that set are frozen, and the same number of positions are admitted from the frozen set, according to the chosen reliability order. The resulting information set $\mathcal I'$ has the same cardinality $K$.

The convolutional precoder is then applied to this **shifted** profile. In particular, a newly frozen position $f$ satisfies $v_f=0$, but its polar input bit need not vanish:

$$u_f=\sum_{\ell=1}^{s}p_\ell v_{f-\ell}\pmod 2.$$

These constraint-dependent input bits are the reason for the reduction in minimum-weight multiplicity described in [5]. The program constructs $\mathcal I'$ first, then uses the same PAC coset calculation described above, with the new frozen positions included in the constraint checks. Thus the PS-PAC run is not an ordinary PAC count followed by an adjustment to the answer: the shifted information set is used throughout the enumeration.

The routine also checks lower row-weight classes introduced by the shift. If their contributions vanish, this does **not** justify moving directly to the next power-of-two distance. The parity-check calculation described in Section 5 is used when necessary to determine whether intermediate weights occur. The PS-PAC $(64,32)$ example below illustrates this point. Ellouze *et al.* [12, 13] also studied minimum-distance and low-weight enumeration for PAC and related pre-transformed polar codes, using polar-coset properties and prefix pruning; their algorithm is not implemented here.

### 4. CRC-polar codes

The CRC-polar calculation follows the minimum-weight codeword formation and CRC conditions developed in [5], using the polar construction of [1]. Let $\mathcal I_{\rm CRC}$ contain the $K$ data positions and the $t$ CRC positions $\mathcal R$. The underlying polar minimum row weight is

$$w_0=\min_{i\in\mathcal I_{\rm CRC}}w(g_i).$$

For each free leading row of weight $w_0$, the program forms $\mathcal K_i$ and enumerates the subsets $\mathcal J$. The balancing set $\mathcal M(\mathcal J)$ gives the required polar input vector with support

$$\mathcal O=\{i\}\cup\mathcal J\cup\mathcal M(\mathcal J).$$

A candidate is rejected if $\mathcal O$ contains a frozen index or if the bits of $\mathcal O$ at the CRC positions differ from the remainder calculated from its data bits. Otherwise it contributes one to the count for the leading coset. The coset counts are added to obtain $A_{w_0}$ for the CRC-polar code. The code uses a direct subset traversal for this test; it does not reproduce every early-pruning step of Algorithm 1 in [5].

If $A_{w_0}>0$, then $d_{\min}=w_0$ and this count is $A_{d_{\min}}$. If $A_{w_0}=0$, the CRC has removed the lowest-weight class, but words of weight $1.5w_0$ or another intermediate weight may remain [5]. The program then needs the separate low-weight search described next.

### 5. Low-weight search using parity-check syndromes

**This part is used only when the structural method does not establish the minimum distance.** It does not run for the usual polar formula or when the minimum-weight class has already been counted exactly.

The search uses a full-rank generator matrix for the selected code, including the PAC precoder or CRC constraints, and derives a parity-check matrix $H$. Let $h_j$ be its $j$th column. An output vector supported on a set $S$ is a codeword exactly when

$$\bigoplus_{j\in S}h_j=0.$$

Divide the $N$ output coordinates into disjoint sets $L$ and $R$. For a subset $T$, write $\sigma(T)=\bigoplus_{j\in T}h_j$. Let $M_L(a,s)$ be the number of subsets of $L$ of size $a$ and syndrome $s$, and define $M_R(b,s)$ in the same way. A weight-$w$ support is a codeword if its two partial syndromes agree. Hence

$$A_w=\sum_{a=\max(0,w-|R|)}^{\min(w,|L|)}\;\sum_s
M_L(a,s)M_R(w-a,s).$$

For each split, the script tabulates the syndrome frequencies on one side and matches them with syndromes from the other side. The frequencies are needed because different subsets can have the same syndrome. Each support has a unique split into $L$ and $R$, so it is counted once. The search tests successive weights until it finds the first nonzero $A_w$. It therefore determines both $d_{\min}$ and $A_{d_{\min}}$, provided all smaller weights have been ruled out.

For example, the $[3,2,2]$ even-parity code has $H=[1\;1\;1]$. Taking $L=\{0\}$ and $R=\{1,2\}$, at weight two the split $(0,2)$ contributes one support and $(1,1)$ contributes two; thus $A_2=3$.

This is an exact **support** enumeration, not a full $2^K$ message enumeration. It is based on parity-check algebra and a syndrome-collision search. Stern's classical low-weight search [11] is relevant to the collision idea, but the present weight-by-weight counter is not Stern's probabilistic algorithm. The particular counting procedure is specified by the expression above rather than attributed to [11].

A search at weight $w$ can still involve many subsets. The following limits apply **to the syndrome search only**:

| Parameter | Default | Meaning |
| --- | ---: | --- |
| `max_exact_weight` | 10 | Highest weight to be tested by the syndrome method |
| `max_syndrome_combos` | 12,000,000 | Maximum subset combinations allowed for either side of a split |
| `max_hash_entries` | 1,000,000 | Maximum number of entries allowed in the syndrome table |

They can be changed in `CONFIG` or on the command line. These values limit computation; they are not target weights supplied by the user. If a required weight or split exceeds a limit, the program reports `NOT CERTIFIED`, rather than claiming a value for $A_{d_{\min}}$. It may also report bounds on $d_{\min}$.

A result contains `d_min`, `A_dmin` and `certified` when called as a Python function. `weight_counts` records weights counted exactly by the support search; `structural_classes` records contributions from particular polar-row weight classes and should **not** be read as a complete weight enumerator. `per_coset=True` gives the individual leading-coset counts for the structural method. The usual screen output prints the same distance and multiplicity information.

## Examples using `CONFIG`

For each example below, replace the `CONFIG` dictionary at the top of **either** script and run it without command-line arguments. The examples use the index and polynomial conventions given above.

### Polar $(128,64)$, crossed RM-polar

```python
CONFIG = dict(
    code_type="polar", N=128, K=64,
    construction="crossed-rm", design_snr_db=3.5,
)
```

Expected: $d_{\min}=16$, $A_{d_{\min}}=94\,488$.

### PAC $(128,64)$, crossed RM-polar

```python
CONFIG = dict(
    code_type="pac", N=128, K=64,
    construction="crossed-rm", design_snr_db=3.5,
    pac_poly="1011011",
)
```

Expected: $d_{\min}=16$, $A_{d_{\min}}=3\,120$ [3, 4].

### PS-PAC $(64,32)$, GA construction

```python
CONFIG = dict(
    code_type="ps-pac", N=64, K=32,
    construction="dega", design_snr_db=2.0,
    alpha=7, pac_poly="1011011111",
)
```

Expected: $A_4=A_5=A_6=A_7=0$ and $d_{\min}=8$, $A_{d_{\min}}=23$ [5]. This example invokes the syndrome search to rule out weights below eight. The polynomial `1011011111` is used in the numerical table of [5]; another passage of that manuscript gives `1011011011`, which does not reproduce the tabulated count.

### CRC-polar $(64,32+5)$

```python
CONFIG = dict(
    code_type="crc-polar", N=64, K=32,
    construction="dega", design_snr_db=2.0,
    crc_poly="x^5+x^3+1",
)
```

Expected: $d_{\min}=8$, $A_{d_{\min}}=64$ [5].

### CRC-polar $(128,64+8)$, PW construction

```python
CONFIG = dict(
    code_type="crc-polar", N=128, K=64,
    construction="pw", crc_poly="0x1D5",
)
```

Expected: $d_{\min}=8$, $A_{d_{\min}}=14$. This example uses the PW profile and the specified CRC polynomial; results with another information set or polynomial need not agree.

### CRC-polar $(32,16+5)$, increased minimum distance

```python
CONFIG = dict(
    code_type="crc-polar", N=32, K=16,
    construction="explicit", info_set="7,11-15,17-31",
    crc_poly="x^5+x^3+1",
)
```

The CRC removes all weight-four words. The syndrome calculation finds $A_4=A_5=0$ and $d_{\min}=6$, $A_{d_{\min}}=14$. This example also shows why the next minimum distance need not be a power of two.

## Execution time

Numba is useful when the coset enumeration involves many combinations. For PAC $(128,64)$ with the crossed RM-polar profile and `pac_poly="1011011"`, an indicative measurement on an x86-64 Linux system was approximately **1.9 s** with pure Python and **0.13 s** with Numba, after compilation. The times depend on the code parameters and processor; these are not Windows measurements.

**The first Numba run can take longer** because Numba compiles the relevant functions into machine code at their first use (just-in-time compilation). Later calls normally reuse compiled code. The script enables caching where supported, although recompilation may be required after changes to the code or input types. The Numba time quoted above excludes the first-call compilation cost.

## References

[1] M. Rowshan, S. H. Dau and E. Viterbo, “On the Formation of Min-Weight Codewords of Polar/PAC Codes and Its Applications,” *IEEE Transactions on Information Theory*, vol. 69, no. 12, pp. 7627–7649, 2023. [doi:10.1109/TIT.2023.3319015](https://doi.org/10.1109/TIT.2023.3319015).

[2] M. Bardet, V.-F. Drăgoi, A. Otmani and J.-P. Tillich, “Algebraic Properties of Polar Codes From a New Polynomial Formalism,” *IEEE International Symposium on Information Theory (ISIT)*, pp. 230–234, 2016. [doi:10.1109/ISIT.2016.7541295](https://doi.org/10.1109/ISIT.2016.7541295).

[3] M. Rowshan and J. Yuan, “Fast Enumeration of Minimum Weight Codewords of PAC Codes,” *IEEE Information Theory Workshop (ITW)*, pp. 255–260, 2022. [doi:10.1109/ITW54588.2022.9965901](https://doi.org/10.1109/ITW54588.2022.9965901).

[4] M. Rowshan and J. Yuan, “On the Minimum Weight Codewords of PAC Codes: The Impact of Pre-Transformation,” *IEEE Journal on Selected Areas in Information Theory*, vol. 4, pp. 487–498, 2023. [doi:10.1109/JSAIT.2023.3312678](https://doi.org/10.1109/JSAIT.2023.3312678).

[5] X. Gu, M. Rowshan and J. Yuan, “CRC-Polar-Inspired PAC Codes: Enumeration and Minimum-Weight Codewords Pruning,” accepted for publication in *IEEE Transactions on Communications*.

[6] M. Rowshan, A. Burg and E. Viterbo, “Polarization-Adjusted Convolutional (PAC) Codes: Sequential Decoding vs List Decoding,” *IEEE Transactions on Vehicular Technology*, vol. 70, no. 2, pp. 1434–1447, 2021. [doi:10.1109/TVT.2021.3052550](https://doi.org/10.1109/TVT.2021.3052550).

[7] P. Trifonov, “Efficient Design and Decoding of Polar Codes,” *IEEE Transactions on Communications*, vol. 60, no. 11, pp. 3221–3227, 2012. [doi:10.1109/TCOMM.2012.081512.110872](https://doi.org/10.1109/TCOMM.2012.081512.110872).

[8] 3GPP TS 38.212, *Multiplexing and Channel Coding*, polar reliability sequence.

[9] Y. Zhou, R. Li, H. Zhang, H. Luo and J. Wang, “Polarization Weight Family Methods for Polar Code Construction,” *IEEE Vehicular Technology Conference (VTC Spring)*, 2018. [doi:10.1109/VTCSpring.2018.8417498](https://doi.org/10.1109/VTCSpring.2018.8417498).

[10] G. He *et al.*, “β-expansion: A Theoretical Framework for Fast and Recursive Construction of Polar Codes,” arXiv:1704.05709, 2017. [arXiv](https://arxiv.org/abs/1704.05709).

[11] J. Stern, “A Method for Finding Codewords of Small Weight,” in *Coding Theory and Applications*, Lecture Notes in Computer Science, vol. 388, pp. 106–113, Springer, 1989. [doi:10.1007/BFb0019850](https://doi.org/10.1007/BFb0019850).

[12] M. Ellouze, R. Tajan, C. Leroux, C. Jégo and C. Poulliat, “Low-Complexity Algorithm for the Minimum Distance Properties of PAC Codes,” *International Symposium on Topics in Coding (ISTC)*, 2023. [doi:10.1109/ISTC57237.2023.10273572](https://doi.org/10.1109/ISTC57237.2023.10273572).

[13] M. Ellouze, R. Tajan, C. Leroux, C. Jégo and C. Poulliat, “Computing the Low-Weight Codewords of Punctured and Shortened Pre-Transformed Polar Codes,” *International Symposium on Topics in Coding (ISTC)*, 2025. [doi:10.1109/ISTC65386.2025.11154614](https://doi.org/10.1109/ISTC65386.2025.11154614).
